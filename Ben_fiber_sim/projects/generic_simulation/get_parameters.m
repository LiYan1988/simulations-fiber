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

%% Plotting parameters
p.fig.spectrum   = subfigure(3, 4, 1, figure(1)); clf;
p.fig.eye        = subfigure(3, 4, 4, figure(2)); clf;
p.fig.const      = subfigure(3, 4, 7, figure(3)); clf;
p.fig.spectrum2  = subfigure(3, 4, 2, figure(4)); clf;
p.fig.eye2       = subfigure(3, 4, 5, figure(5)); clf;
p.fig.const2     = subfigure(3, 4, 8, figure(6)); clf;

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

% Free parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data rate etc.
p.f_symb        = 10e9; % The symbol rate [Baud]
p.samp_per_symb = 32;   % Samples per symbol
p.N_symb        = 256;  % Number of symbols

%% Input signal
%p.modulation.format = 'OOK'; p.bits_per_symbol = 1;
%p.modulation.format = 'BPSK'; p.bits_per_symbol = 1;
%p.modulation.format = 'QPSK'; p.bits_per_symbol = 2;
p.modulation.format = '16-QAM'; p.bits_per_symbol = 4;
p.data = round(rand(p.bits_per_symbol, p.N_symb));

p.light.power_dBm   = 15; % Mean power [dBm]

e_levels           = [0 acos(1/3)/pi acos(-1/3)/pi 1]; % Modulation levels
%e_filter.type      = 'none';     % No filtering applied
%e_filter.type      = 'ideal_bp';     % Ideal (rectangular) filter
e_filter.type      = 'gaussian';     % A Gaussian filter
e_filter.bandwidth = 0.65*p.f_symb;  % Filter bandwidth [Hz]
e_filter.bw_spec   = 'HWHM';         % Half-width, half maximum

%p.electric1.data   = p.data(1:1, :); % Data for electric source
p.electric1.data   = p.data(1:2, :); % Data for electric source
p.electric1.levels = e_levels;
p.electric1.filter = e_filter;

%p.electric2.data   = p.data(2:2, :); % Data for electric source
p.electric2.data   = p.data(3:4, :); % Data for electric source
p.electric2.levels = e_levels;
p.electric2.filter = e_filter;

%% Transmission parameters
% An SMF
% Transmission distance [m]
p.smf.L = 20e3; % Transmission distance [m]
% Group-velocity dispersion [s^2/m]
p.smf.beta2 = 17*p.conv.D_to_beta2;
% Attenuation [1/m]
p.smf.alpha = 0.25*p.conv.att_DB_to_att;
% Kerr nonlinearity [1/W/m]
p.smf.gamma = 1.3e-3;
% Step size for NLSE solver
p.steps_per_L_NL = 200;

% Derived parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
