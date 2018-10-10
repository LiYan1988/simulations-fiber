function p = get_parameters(p)
% Set up parameters for a generic simulation.
%
% INPUTS:
%    p : The parameter struct from main
% OUTPUTS:
%    p : The parameter struct
%
% Pontus Johannisson, 2010-04-19
% This software is distributed under the terms of the GNU General
% Public License version 2

%% Flags
p.flag.time_seed = 0; % Flag for using current time to seed RNG (1 = use current time; 0 = use fixed seed)
p.flag.dual_pol  = 0; % Flag for dual-polarization operation
p.flag.gpu       = 0; % Flag to use GPU acceleration
p.flag.edc_es    = 1; % Flag to perform edc at every span.

p.flag.osa       = 1; % Flag for plotting spectra
p.flag.rx_perf   = 1; % Flag for receiver performance monitoring

p.flag.b2b       = 0; % Flag for back-to-back operation
p.flag.psa       = 0; % Flag for PSA operation


%% Constants
p.const.c      = 2.99792458e8;             % Speed of light [m/s]
p.const.q      = 1.6021773e-19;            % Electron charge [C]
p.const.h      = 6.626076e-34;             % Planck's constant [Js]
p.const.lambda = 1.55e-6;                  % Wavelength [m]
p.const.nu     = p.const.c/p.const.lambda; % Frequency [Hz]

%% Conversions
p.conv.D_to_beta2    = -p.const.lambda^2/2/pi/p.const.c*1e-6; % * D [ps/nm/km]
p.conv.att_DB_to_att = log(10)/10/1e3; % * alpha_dB [dB/km]
p.conv.nm_to_Hz      = p.const.c/p.const.lambda^2*1e-9; % * BW [nm]

%% Free parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Random number generation
if p.flag.time_seed
    if ~isfield(p, 'timestamp')
        p.timestamp = now();
    end
    p.rng.seed = mod(p.timestamp*1e6, 2^32-1); % The time vector gives number of whole days since epoch. RNG seed takes 32-bit integers. Multiplying by 1e6 means that the rng seed changes every ~86ms
else
    p.rng.seed = 1;
end
rng(p.rng.seed); 

%% Log file for results
% Path to save results to
% p.file.path = '/c3se/NOBACKUP/users/foob/nlse_simulations/results/'; % Path to save files to
p.file.path = '../../results/'; % For local testing 
if ~exist(p.file.path, 'dir')
    mkdir(p.file.path)
end

log_id = getCurrentJob();
if ~isempty(log_id)
    log_id = ['_', num2str(log_id.ID)];
else
    log_id = ''; % In case this is not run via batch processing
end
[func_stack, ~] = dbstack();
p.log.filename = sprintf('log_%s%s.csv', func_stack(end).file(1:end-2), log_id); % Log file

% Open log file
if isfield(p.log, 'filename')
    p.log.fid = fopen(p.log.filename, 'a'); % Open Log file in append mode
    if p.log.fid == -1
        warning('Could not open log file')
    end
end

% Saving output of every Xth span
p.file.save_span_incr = 5;

%% Data rate and Link parameters
%p.wdm_power= [-12 -10 -8 -6 -4 -2 0 2 4 6 8 10] ; %Defines the total power of the WDM signal, When off, the defined power per channel is applied.

p.samp_per_symb = 4;   % Samples per symbol
p.N_symb        = 2^14;  % Number of bits

p.f_symb        = [10e9]; % The symbol rate [Baud]
p.link.modulations      = [1];%Set the modulation format for each channel. 1 OOK, 2 BPSK, 3 QAM, 4 16QAM
p.power_array        = [-20:2:20]; %Define the power difference between WDM channels.

% p.f_symb        = [10e9 10e9 10e9 10e9 10e9 30e9 10e9 10e9 10e9 10e9 10e9]; % The symbol rate [Baud] 
% p.link.modulations      = [1 1 1 1 1 4 1 1 1 1 1]; %1 OOK, 2 BPSK, 3 QAM, 4 16QAM
% p.power_array        = [0 0 0 0 0 0 0 0 0 0 0]; %Define the power difference between WDM channels.

p.link.N_span        = 5;    % Number of spans in the link
p.link.N_chan        = length(p.f_symb);    % Number of WDM channels
p.link.edfa_nf_db    = 4;    % EDFA noise figure [dB]


%% Transmitter paramters
%p.tx.filter_bandwidth_hz = p.f_symb*2;
p.tx.chan_spacing_hz = 50e9; % Space between WDM channels [Hz]
p.tx.linewidth = 00e3; % Laser linewidth [Hz]
p.tx.snr_db = 2000; % Input noise, this is not exatly SNR [dB]. A model for transmitter imperfections.

p.w_index.sig = zeros(1, p.link.N_chan);
for chan_ind = 1:p.link.N_chan
    p.w_index.sig(chan_ind) = 2*pi*((2*chan_ind-1-p.link.N_chan)/2)*p.tx.chan_spacing_hz; %For Multiplexing the WDM signal
end

%% Transmission fiber parameters
% An SMF
% Span length [m]
p.smf.L = 82e3;
% Dispersion parameter [s/m^2]
p.smf.D = 17e-6;
% Group-velocity dispersion [s^2/m]
p.smf.beta2 = p.smf.D*1e6*p.conv.D_to_beta2;
% Attenuation [dB/m]
p.smf.att_db = 0.2e-3;
% Attenuation [1/m]
p.smf.alpha = p.smf.att_db*1e3*p.conv.att_DB_to_att;
% Kerr nonlinearity [1/W/m]
p.smf.gamma = 1.4e-3;
% Step size for NLSE solver
p.steps_per_L_NL = 100;

%% Dispersion management parameters (important for NLC with PSA)
p.psa.pre_comp_L = 0e3; % Dispersion pre-compensation length [m]

%% Receiver parameters

p.rx.filter_bandwidth_hz = p.f_symb*2; %RX filter Bandwidth
p.rx.edc_L = p.smf.L*p.link.N_span; %Total Electronic Dispertion Compensation. ONLY WORKS if p.flag.edc_es is 0.
p.rx.responsivity  = 1; % Photodiode responsivity [A/W] (Value between 0 and 1)
p.rx.thermal_noise = 1.81989010657e-11; % Single-sided thermal noise spectral density in photodiode [A/sqrt(Hz)] - Assumes operation at 300 K with 50-ohm load resistor
p.rx.linewidth = 00e3; % Laser linewidth [Hz]

% Equalizer properties

p.rx.vv_block_length = 20; 
p.rx.os = 1; 
p.rx.equalizer.NTap             = 25; % Number of equalizer taps
p.rx.equalizer.StepSize_CMA     = 1e-3; % Step size of the CMA pre-convergence equalizer
p.rx.equalizer.RefPower_CMA     = 1.32; % CMA reference power
p.rx.equalizer.Iterations_CMA   = 4; % Number of CMA iterations
p.rx.equalizer.Iterations_RDCMA = 0; % Number of Radial-decision CMA iterations
p.rx.equalizer.SamplesPerSymbol = p.rx.os; % Number of samples per symbol
p.rx.equalizer.SymbolRate       = p.f_symb; % Symbol rate 
p.rx.equalizer.BlockSize        = 40; % Phase estimator block size
p.rx.equalizer.fftSize_FOE      = 2^8; % FFT size for frequency offset estimation
p.rx.equalizer.TestAngles_DDLMS = 128; % Number of test angles for phase estimation



%% Derived parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.f_samp   = max(p.f_symb)*p.samp_per_symb; % Sampling rate [Hz]
p.N_samp   = p.N_symb*p.samp_per_symb;
p.dt_symb  = 1/max(p.f_symb);               % Symbol interval [s]
p.dt_samp  = 1/p.f_samp;               % Sampling interval [s]
p.T_period = p.N_samp*p.dt_samp;
p.t        = 0:p.dt_samp:(p.N_samp - 1)*p.dt_samp;
p.f        =   1 /p.T_period*[(0:p.N_samp/2 - 1) (-p.N_samp/2:-1)];
p.omega    = 2*pi/p.T_period*[(0:p.N_samp/2 - 1) (-p.N_samp/2:-1)];

% Set up indices in the middle of the bit slot
p.idx_symb = (p.samp_per_symb/2 + 1):p.samp_per_symb:p.N_samp;

