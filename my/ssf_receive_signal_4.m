clc;
clear;
close all;

% Wrap chromatic dispersion compensation. 
% Try compensate for constellation rotation due to nonlinearity.

%% Fiber Parameters
% -------------- Primary parameters
param.fmax = 2*pi*400*1e9; % [Hz]
param.fn = 2^22; % number of spectrum points

param.span_length = 100; % [km], span length
param.beta2 = -2.1683e-23; % [s^2/km], GVD, D=17 [ps/ns/km]
param.gamma = 1.27; % [(W*km)^-1], nonlinear coefficient of SMF
param.alpha = 0*log(10)*0.2/10; % [1/km] in linear, 0.2 dB/km, positive number
param.zn = 1000; % number of steps per span

%% Channel Parameters
% Channel specific parameters, n channels should have n sets of parameters
N = 5; % number of channels, should be an odd number
param.bandwidth_channel = 10*1e9*ones(1, N); % [Hz]
param.bandwidth_channel((N-1)/2+1) = 32*1e9; % CUT
param.constellation_size = 2*ones(1, N); % 2=BPSK; 4=QPSK; 8=8QAM; 16=16QAM; 32=32QAM; 64=64QAM
% param.constellation_size = [2, 4, 8, 16, 32];
param.constellation_size((N-1)/2+1) = 16;
param.spectrum_grid_size = 50*1e9; % [Hz], spectrum grid size
param.center_frequency_channel = param.spectrum_grid_size*(linspace(0, N-1, N)-(N-1)/2);
param.power_channel_time = 10^(-6/10)/1e3*ones(N, 1); % [W], power of channel in time domain, in contrast to the frequency domain PSD measured in W/Hz

% Filter parameters of each channel, assume raised cosine filters
param.roll_off_filter = 0.1*ones(1, N); % Roll-off factor of raised cosine filter
param.symbol_in_filter = 100*ones(1, N); % length of impulse response in symbol
param.shape_filter = cell(1, N); % shape of raised
param.shape_filter(:) = {'normal'};

% Random seed
param.random_seed = 2394759; % input to rng

%% Generate Signal
param = generate_signals(param);

%% Split Step Fourier
param = split_step_single_polarization(param);
% plot_current_signal(param, 'linear')

%% Plot Transmitted Signal Before and After Chromatic Dispersion Compensation
% cidx = (N-1)/2+1; % index of the channel to display
cidx = 3;
[xt_dc, xf_dc, xt] = dispersion_compensation(param, cidx);
scatterplot(xt_dc(param.delay_filter_channel(cidx)+1:end-param.delay_filter_channel(cidx)), ...
    param.sample_per_symbol(cidx), 0)

%%