clc;
clear;
close all;

%% Fiber Parameters
% -------------- Primary parameters
param.fmax = 2*pi*400*1e9; % [Hz]
param.fn = 2^22; % number of spectrum points

param.span_length = 10; % [km], span length
param.beta2 = -2.1683e-05; % [ns^2/km], GVD
param.gamma = 1.27; % [(W*km)^-1], nonlinear coefficient of SMF
param.alpha = log(10)*0.2/10; % [1/km] in linear, 0.2 dB/km, positive number
param.zn = 100; % number of steps per span

%% Channel Parameters
% % Channel specific parameters, n channels should have n sets of parameters
% param.bandwidth_channel = [10*1e9, 10*1e9, 32*1e9, 10*1e9, 10*1e9]; % [Hz]
% param.constellation_size = [2, 2, 16, 2, 2]; % 2=BPSK; 4=QPSK; 8=QAM; 16=16QAM; 32=32QAM; 64=64QAM
% param.center_frequency_channel = [-100*1e9, -50*1e9, 0*1e9, 50*1e9, 100*1e9];
%
% % Filter parameters of each channel, assume raised cosine filters
% param.roll_off_filter = [0.1, 0.1, 0.1, 0.1, 0.1]; % Roll-off factor of raised cosine filter
% param.symbol_in_filter = [10, 10, 10, 10, 10]; % length of impulse response in symbol
% param.shape_filter = {'normal', 'normal', 'normal', 'normal', 'normal'}; % shape of raised

% Channel specific parameters, n channels should have n sets of parameters
N = 3; % number of channels, should be an odd number
param.bandwidth_channel = 10*1e9*ones(1, N); % [Hz]
% param.bandwidth_channel((N-1)/2+1) = 32*1e9; % CUT
param.constellation_size = 2*ones(1, N); % 2=BPSK; 4=QPSK; 8=QAM; 16=16QAM; 32=32QAM; 64=64QAM
% param.constellation_size((N-1)/2+1) = 16;
param.center_frequency_channel = 50*1e9*(linspace(0, N-1, N)-(N-1)/2);
param.power_channel_time = 10^(-6/10)/1e3*ones(N, 1); % [W], power of channel in time domain, in contrast to the frequency domain PSD measured in W/Hz

% Filter parameters of each channel, assume raised cosine filters
param.roll_off_filter = 0.1*ones(1, N); % Roll-off factor of raised cosine filter
param.symbol_in_filter = 10*ones(1, N); % length of impulse response in symbol
param.shape_filter = cell(1, N); % shape of raised
param.shape_filter(:) = {'normal'};

% Random seed
param.random_seed = 2394759; % input to rng

%% Generate Signal
param = generate_signals(param);

%% Plot signal
plot_current_signal(param, 'linear');

P_t = norm(param.data_mod_t_in)^2/param.fn;
P_f = norm(param.data_mod_f_in)^2/(param.fn*param.dt);

%% Split Step Fourier
param = split_step_single_polarization(param);
plot_current_signal(param, 'linear')

%% Receive Signal
x = param.data_mod_t_current;

%%
% f2 = (param.f_plot<25)&(param.f_plot>-25);
% xf2 = f2.*param.data_mod_f_current;
%
% xt2 = ift(xf2, param.fn, param.df);
%
% figure;
% subplot(2, 1, 1)
% plot(param.t_plot, abs(xt2).^2)
% subplot(2, 1, 2)
% plot(param.f_plot, abs(xf2).^2)
%
% %%
% filt2 = param.filter_tx_channel{2};
% xr2 = upfirdn(xt2, filt2, 1, 1);
%
% plot(abs(xr2).^2)

%% Back to Back Receiver of BPSK signal
clc;
close all;
% Transmitted waveform
xt2 = param.data_mod_t_channel{2};
xt2 = xt2/norm(xt2);
% plot(param.t_plot, xt2)

% square-root, raised cosine filter used in both transmitter and receiver
filt2 = param.filter_tx_channel{2};

% pass signal through the filter
xr2 = upfirdn(xt2, filt2, 1, 1);
xr2 = xr2/norm(xr2);
% length difference between input and output of the filter is the length of
% filter minus 1
length_difference = length(xr2)-length(xt2);

% The input and output singals are very similar to each other. The output
% signal is a shifted and padded version of the input. Use
% cross-correlation to calculate the time shift between them.
% xcorr calculates the similarity between xt2 and shifted xr2
[r, lags] = xcorr(xt2, xr2);
% plot(lags, r)
[~, idx] = max(r);
lag_max = lags(idx);
% lag_max is half the length of the filter, (length(filt2)-1)/2
xr2 = circshift(xr2, lag_max);
xr2 = xr2(1:end-length_difference);

% plot input and output of the filter together to align them in time
% figure;
% hold on;
% plot(xt2)
% plot(xr2)

% Now xr2 and xt2 are the same length. Down sample xr2 and plot
% constellation.
pad_length = param.sample_per_symbol(2)*ceil(length(xr2)/param.sample_per_symbol(2))-length(xr2);
xr2 = [xr2; zeros(pad_length, 1)];
xr2 = reshape(xr2, param.sample_per_symbol(2), []);
% p = xr2>0;
% xr2p = xr2.*p;
% plot(mean(xr2p, 2))
[~, idx] = max(mean(xr2.*(xr2>0), 2));
xr2d = xr2(idx, :)';
xr2d = (xr2d-mean(xr2d))/(max(xr2d)-min(xr2d));

% figure;
% hold on;
% plot(xr2d)
% plot(xt2d-0.5)

scatterplot(xr2d)
hold on
scatterplot(xt2d-0.5)

