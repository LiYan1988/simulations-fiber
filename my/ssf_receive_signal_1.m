clc;
clear;
close all;

% Receive signal in back-to-back scenario

%% Fiber Parameters
% -------------- Primary parameters
param.fmax = 2*pi*400*1e9; % [Hz]
param.fn = 2^20; % number of spectrum points

param.span_length = 10; % [km], span length
param.beta2 = -2.1683e-05; % [ns^2/km], GVD
param.gamma = 1.27; % [(W*km)^-1], nonlinear coefficient of SMF
param.alpha = log(10)*0.2/10; % [1/km] in linear, 0.2 dB/km, positive number
param.zn = 100; % number of steps per span

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
param.power_channel_time = [1:N]*1e-3; % [W], power of channel in time domain, in contrast to the frequency domain PSD measured in W/Hz

% Filter parameters of each channel, assume raised cosine filters
param.roll_off_filter = 0.1*ones(1, N); % Roll-off factor of raised cosine filter
param.symbol_in_filter = 100*ones(1, N); % length of impulse response in symbol
param.shape_filter = cell(1, N); % shape of raised
param.shape_filter(:) = {'normal'};

% Random seed
param.random_seed = 2394759; % input to rng

%% Generate Signal
param = generate_signals(param);

%% Plot signal
% plot_current_signal(param, 'linear');

P_t = norm(param.data_mod_t_in)^2/param.fn;
P_f = norm(param.data_mod_f_in)^2/(param.fn*param.dt);

%% Split Step Fourier
param = split_step_single_polarization(param);
% plot_current_signal(param, 'linear')

%% Plot Transmitted Signal
cidx = 5; % index of the channel to display
% create a mask in at the baseband
% param.mask_f = (param.f_plot<25)&(param.f_plot>-25);
% downconvert selected channel and low-pass filter it

% plot_current_signal(param, 'linear')

%%
x = param.data_mod_t_in.*exp(1i*2*pi*param.center_frequency_channel(cidx).*param.t);
xf = ft(x, param.df);
xf = xf.*param.f_mask;
% plot(param.f_plot, abs(xf).^2)
xt = ift(xf, param.df);

%%
s = 1;
delay_filter = param.symbol_in_filter(cidx)/2*param.sample_per_symbol(cidx);
padding_length = ceil(length(xt)/param.sample_per_symbol(cidx))*...
    param.sample_per_symbol(cidx)-length(xt);
xp = reshape([xt; zeros(padding_length, 1)], ...
    param.sample_per_symbol(cidx), []);
scatter(real(xp(s, delay_filter+1:end-delay_filter)), imag(xp(s, delay_filter+1:end-delay_filter)))

%%
% x = param.data_mod_t_channel{cidx};
% % plot(param.f_plot, abs(ft(x, param.df)).^2)
% y = x.*exp(1i*2*pi*param.center_frequency_channel(cidx).*param.t);
% yf = ift(y, param.df);
% plot(param.f_plot, abs(yf).^2)