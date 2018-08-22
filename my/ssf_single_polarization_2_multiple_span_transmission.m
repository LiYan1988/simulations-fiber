clc;
clear;
close all;

% Add third order dispersion

%% Fiber Parameters
% -------------- Primary parameters
param.fmax = 2*pi*200*1e9; % [Hz]
param.fn = 2^17; % number of spectrum points

param.span_length = 82; % [km], span length
param.beta2 = -2.1683e-23; % [s^2/km], GVD, D=17 [ps/ns/km]
S = 0.06*1e6; % [s/(m^2*km)], third order dispersion
param.wavelength = 1550*1e-9; % [m], reference wavelength
param.light_speed = 2*1e8; % [m/s], speed of light in fiber
param.beta3 = (S - 4*pi*param.light_speed/param.wavelength^3*param.beta2)*...
    (param.wavelength^2/(2*pi*param.light_speed))^2; % [s^3/km], third order dispersion
% d3=1i*beta3/6*(2*pi*FF).^3;
param.gamma = 1.27; % [(W*km)^-1], nonlinear coefficient of SMF
param.alpha = log(10)*0.2/10; % [1/km] in linear, 0.2 dB/km, positive number
param.zn = 1000; % number of steps per span

%% FBG
fbg_beta2 = -param.beta2;
fbg_beta3 = -param.beta3;
fbg1_length = 80; % [km]
fbg1_dispersion = exp(0.5*1i*fbg_beta2*param.f.^2*fbg1_length).*...
    exp(1i/6*fbg_beta3*param.f.^3*fbg1_length);

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
param.power_channel_time = 10^(-3/10)/1e3*ones(N, 1); % [W], power of channel in time domain, in contrast to the frequency domain PSD measured in W/Hz

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
cidx = (N-1)/2+1; % index of the channel to display
% cidx = 2;

[xt_dc, ~, ~] = dispersion_compensation(param.data_mod_t_current, param, cidx);
scatterplot(xt_dc(param.delay_filter_channel(cidx)+1:end-param.delay_filter_channel(cidx)), ...
    param.sample_per_symbol(cidx), param.shift_channel_time(cidx))

%% Compensate nonlinear rotation
% This can only compensate for a very little part of the constellation
% rotation because most of the rotation is due to XPM
% Moreover, a simple constellation rotation cannot compensate for the
% nonlinear phase noise. Backpropagation is needed.

% if param.alpha>0
%     param.span_length_effective = (1-exp(-param.alpha*param.span_length))/param.alpha;
% elseif param.alpha==0
%     param.span_length_effective = param.span_length;
% end
% 
% nonlinear_phase = -1i*param.gamma*param.span_length_effective*exp(param.alpha*param.span_length);
% 
% % nonlinear_phase = -param.hhz*param.zn;
% 
% xt_nc = xt_dc.*exp(nonlinear_phase.*abs(xt_dc).^2);
% scatterplot(xt_nc(param.delay_filter_channel(cidx)+1:end-param.delay_filter_channel(cidx)), ...
%     param.sample_per_symbol(cidx), param.shift_channel_time(cidx))

%% Backpropagation
param.data_mod_t_in_bp = param.data_mod_t_current;
param = back_propagation_split_step_single_polarization(param);

%%
cidx = (N-1)/2+1; % index of the channel to display
cidx = 3;

[xt_dc, xf_dc, xt] = dispersion_compensation(param.data_mod_t_current_bp, param, cidx);
scatterplot(xt(param.delay_filter_channel(cidx)+1:end-param.delay_filter_channel(cidx)), ...
    param.sample_per_symbol(cidx), param.shift_channel_time(cidx))
