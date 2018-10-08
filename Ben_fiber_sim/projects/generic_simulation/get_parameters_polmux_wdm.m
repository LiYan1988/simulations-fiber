function p = get_parameters_polmux_wdm(p)
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

%% Plotting parameters
% p.fig.spectrum   = subfigure(3, 4, 1, figure(1)); clf;
% p.fig.eye        = subfigure(3, 4, 4, figure(2)); clf;
% p.fig.const      = subfigure(3, 4, 7, figure(3)); clf;
% p.fig.spectrum2  = subfigure(3, 4, 2, figure(4)); clf;
% p.fig.eye2       = subfigure(3, 4, 5, figure(5)); clf;
% p.fig.const2     = subfigure(3, 4, 8, figure(6)); clf;

%% Flags
p.flag.time_seed = 0; % Flag for using current time to seed RNG (1 = use current time; 0 = use fixed seed)
p.flag.dual_pol  = 0; % Flag for dual-polarization operation
p.flag.gpu       = 1; % Flag to use GPU acceleration

p.flag.osa       = 0; % Flag for plotting spectra
p.flag.rx_perf   = 1; % Flag for receiver performance monitoring

p.flag.b2b       = 0; % Flag for back-to-back operation

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
log_id = getCurrentJob();
if ~isempty(log_id)
    log_id = ['_', num2str(log_id.ID)];
else
    log_id = ''; % In case this is not run via batch processing
end
[func_stack, ~] = dbstack();
p.log.filename = ['../results/log_', func_stack(end).file(1:end-2), log_id, '.csv']; % Log file

% Open log file
if isfield(p.log, 'filename')
    p.log.fid = fopen(p.log.filename, 'a'); % Open Log file in append mode
    if p.log.fid == -1
        warning('Could not open log file')
    end
end

%% Data rate etc.
p.f_symb        = 10e9; % The symbol rate [Baud]
p.samp_per_symb = 128;   % Samples per symbol
p.N_symb        = 2^12;  % Number of symbols

%% Link parameters
p.link.N_chan        = 1;    % Number of WDM channels
p.link.pwr_per_chan  = -8;   % Power per channel [dBm]
p.link.N_span        = 1;    % Number of spans in the link
p.link.edfa_nf_db    = 5;    % EDFA noise figure [dB]

%% Input signal
% p.modulation.format = 'OOK'; p.bits_per_symbol = 1;
% p.modulation.format = 'BPSK'; p.bits_per_symbol = 1;
p.modulation.format = 'QPSK'; p.bits_per_symbol = 2;
% p.modulation.format = '16-QAM'; p.bits_per_symbol = 4;
p.data = round(rand(p.bits_per_symbol, p.N_symb, (p.flag.dual_pol>0)+1, p.link.N_chan));

p.light.power_dBm   = 15; % Mean power [dBm]

if strcmpi(p.modulation.format , '16-qam')
    e_levels       = [0 acos(1/3)/pi acos(-1/3)/pi 1]; % Modulation levels
elseif strcmpi(p.modulation.format, 'qpsk')
    e_levels       = [0, 1];
% elseif strcmpi(p.modulation.format, 'bpsk')
%     e_levels       = [0, 1, 0.5];
% elseif strcmpi(p.modulation.format, 'ook')
%     e_levels       = [0.5, 1];
else
    error('Could not find modulation format')
end

%e_filter.type      = 'none';     % No filtering applied
%e_filter.type      = 'ideal_bp';     % Ideal (rectangular) filter
e_filter.type      = 'gaussian';     % A Gaussian filter
e_filter.bandwidth = 0.65*p.f_symb;  % Filter bandwidth [Hz]
e_filter.bw_spec   = 'HWHM';         % Half-width, half maximum

p.electric.levels = e_levels;
p.electric.filter = e_filter;

%% Transmitter paramters
p.tx.filter_bandwidth_hz = p.f_symb*2; % Bandwidth of the transmitter MUX filter [Hz]
p.tx.chan_spacing_hz = 50e9; % Space between WDM channels [Hz]

p.w_index.sig = zeros(1, p.link.N_chan);
for chan_ind = 1:p.link.N_chan
    p.w_index.sig(chan_ind) = 2*pi*((2*chan_ind-1-p.link.N_chan)/2)*p.tx.chan_spacing_hz;
end

%% Transmission fiber parameters
% An SMF
% Transmission distance [m]
p.smf.L = 100e3; % Transmission distance [m]
% Dispersion parameter [s/m^2]
p.smf.D = 17e-6;
% Group-velocity dispersion [s^2/m]
p.smf.beta2 = p.smf.D*1e6*p.conv.D_to_beta2;
% Attenuation [dB/m]
p.smf.att_db = 0.25e-3;
% Attenuation [1/m]
p.smf.alpha = p.smf.att_db*1e3*p.conv.att_DB_to_att;
% Kerr nonlinearity [1/W/m]
p.smf.gamma = 1.3e-3;
% Step size for NLSE solver
p.steps_per_L_NL = 200;

%% Receiver parameters
p.rx.filter_bandwidth_hz = p.f_symb*2; % Bsndwidth of the receiver De-MUX filter [Hz]
p.rx.edc_L = p.smf.L*p.link.N_span; % Length of fiber to perform electric dispersion compensation for [m]
p.rx.responsivity  = 1; % Photodiode responsivity [A/W] (Value between 0 and 1)
p.rx.thermal_noise = 1.81989010657e-11; % Single-sided thermal noise spectral density in photodiode [A/sqrt(Hz)] - Assumes operation at 300 K with 50-ohm load resistor

p.rx.vv_block_length = 64;
p.rx.os = 2; % Receiver DSP oversampling

p.rx.equalizer.NTap             = 17; % Number of equalizer taps
p.rx.equalizer.StepSize_CMA     = 1e-3; % Step size of the CMA pre-convergence equalizer
p.rx.equalizer.RefPower_CMA     = 1; % CMA reference power
p.rx.equalizer.Iterations_CMA   = 4; % Number of CMA iterations
p.rx.equalizer.Iterations_RDCMA = 0; % Number of Radial-decision CMA iterations
p.rx.equalizer.StepSize_DDLMS   = 1e-6; % Step size for DD-LMS equalizer
p.rx.equalizer.Iterations_DDLMS = 2; % Number of DD-LMS iterations
p.rx.equalizer.SamplesPerSymbol = p.rx.os; % Number of sampler per symbol
p.rx.equalizer.SymbolRate       = p.f_symb; % Symbol rate 
p.rx.equalizer.BlockSize        = 65; % Phase estimator block size
p.rx.equalizer.fftSize_FOE      = 2^8; % FFT size for frequency offset estimation
p.rx.equalizer.TestAngles_DDLMS = 128; % Number of test angles for phase estimation
p.rx.equalizer.Constellation = [-1-1j 1-1j 1+1j -1+1j];
p.rx.equalizer.Constellation = single(p.rx.equalizer.Constellation/sqrt(mean(abs(p.rx.equalizer.Constellation).^2)));

%% Derived parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.f_samp   = p.f_symb*p.samp_per_symb; % Sampling rate [Hz]
p.N_samp   = p.N_symb*p.samp_per_symb;
p.dt_symb  = 1/p.f_symb;               % Symbol interval [s]
p.dt_samp  = 1/p.f_samp;               % Sampling interval [s]
p.T_period = p.N_samp*p.dt_samp;
p.t        = 0:p.dt_samp:(p.N_samp - 1)*p.dt_samp;
p.f        =   1 /p.T_period*[(0:p.N_samp/2 - 1) (-p.N_samp/2:-1)];
p.omega    = 2*pi/p.T_period*[(0:p.N_samp/2 - 1) (-p.N_samp/2:-1)];

% Set up indices in the middle of the bit slot
p.idx_symb = (p.samp_per_symb/2 + 1):p.samp_per_symb:p.N_samp;
