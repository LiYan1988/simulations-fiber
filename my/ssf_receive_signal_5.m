clc;
clear;
close all;

% Square-root RRC filter at transmitter and receiver

%% Fiber Parameters
% -------------- Primary parameters
param.fmax = 2*pi*500*1e9; % [Hz]
param.fn = 2^19; % number of spectrum points

param.span_length = 82; % [km], span length
param.beta2 = -2.1683e-23; % [s^2/km], GVD, D=17 [ps/ns/km]
S = 0.06*1e6; % [s/(m^2*km)], third order dispersion
param.wavelength = 1550*1e-9; % [m], reference wavelength
param.light_speed = 2*1e8; % [m/s], speed of light in fiber
param.beta3 = (S-4*pi*param.light_speed/param.wavelength^3*param.beta2)*...
    (param.wavelength^2/(2*pi*param.light_speed))^2; % [s^3/km], third order dispersion
% d3=1i*beta3/6*(2*pi*FF).^3;
param.gamma = 1.27; % [(W*km)^-1], nonlinear coefficient of SMF
param.alpha = log(10)*0.2/10; % [1/km] in linear, 0.2 dB/km, positive number
param.dz = 0.1; % [km]

%% Channel Parameters
% Channel specific parameters, n channels should have n sets of parameters
N = 11; % number of channels, should be an odd number
param.bandwidth_channel = 10*1e9*ones(1, N); % [Hz]
param.bandwidth_channel((N-1)/2+1) = 32*1e9; % CUT
param.constellation_size = 2*ones(1, N); % 2=BPSK; 4=QPSK; 8=8QAM; 16=16QAM; 32=32QAM; 64=64QAM
% param.constellation_size = [2, 4, 8, 16, 32];
param.constellation_size((N-1)/2+1) = 16;
param.spectrum_grid_size = 50*1e9; % [Hz], spectrum grid size
param.center_frequency_channel = param.spectrum_grid_size*(linspace(0, N-1, N)-(N-1)/2);
param.power_channel_time = 10^(-3/10)/1e3*ones(N, 1); % [W], power of channel in time domain, in contrast to the frequency domain PSD measured in W/Hz

% Filter parameters of each channel, assume raised cosine filters
param.roll_off_filter = 0.2*ones(1, N); % Roll-off factor of raised cosine filter
param.symbol_in_filter = 100*ones(1, N); % length of impulse response in symbol
param.shape_filter = cell(1, N); % shape of raised
param.shape_filter(:) = {'sqrt'};

% Random seed
param.random_seed = 2394759; % input to rng

%% Generate Signal
param = generate_signals(param);

%% Propagation through a link
% 4 spans of 82 km SSMF followed by DCF with 80 km SSMF-equivalent
% compensation, plus one span of 82 km SSFM followed by DCF with 90 km
% SSMF-equivalent compensation

% FBG parameters
fbg_beta2 = -param.beta2;
fbg_beta3 = -param.beta3;
fbg_length_1 = 80; % [km]
fbg_length_2 = 90; % [km]

for n=1:4
    % Split Step Fourier in SSMF
    param = split_step_single_polarization(param);
    
    % FBG
    [signal_t_out, signal_f_out] = fbg_propagation(...
        param.data_mod_t_current, fbg_beta2, fbg_beta3, fbg_length_1, param);
    
    param.data_mod_t_current = signal_t_out;
    param.data_mod_f_current = signal_f_out;
end

% one span of 82 km SSFM, followed by DCF with 90 km SSMF-equivalent
% compression
param = split_step_single_polarization(param);
[signal_t_out, signal_f_out] = fbg_propagation(...
    param.data_mod_t_current, fbg_beta2, fbg_beta3, fbg_length_2, param);
param.data_mod_t_current = signal_t_out;
param.data_mod_f_current = signal_f_out;

%% Signal after FBG
cidx = (N-1)/2+1; % index of the channel to display

% Plot all the channels
for cidx=(N-1)/2+1
    % compensate for residual dispersion
    [xt_dc, ~, ~] = dispersion_compensation(param.data_mod_t_current, cidx,...
        param.beta2, param.beta3, 0, param);
    scatterplot(...
        xt_dc(param.delay_filter_channel(cidx)+1:...
        end-param.delay_filter_channel(cidx)), ...
        param.sample_per_symbol(cidx), param.shift_channel_time(cidx))
    title(sprintf('Channel %d', cidx))
end

%% Signal clusters
close all;
clc;
signal = zeros(size(xt_dc(param.delay_filter_channel(cidx)+1:...
        end-param.delay_filter_channel(cidx)), 1), 2);
signal(:, 1) = real(xt_dc(param.delay_filter_channel(cidx)+1:...
        end-param.delay_filter_channel(cidx)));
signal(:, 2) = imag(xt_dc(param.delay_filter_channel(cidx)+1:...
        end-param.delay_filter_channel(cidx)));

signal = downsample(signal, param.sample_per_symbol(cidx), ...
    param.shift_channel_time(cidx));

% Find center of points 
opts = statset('UseParallel', true);
[idx,C,sumd,D] = kmeans(signal, 16, 'Display', 'final', ...
    'maxiter', 1000, 'Replicates', 64, 'Options', opts);

figure;
hold on;
box on;
grid on;
plot(signal(:, 1), signal(:, 2), '.')
plot(C(:, 1), C(:, 2), 'x')

centers = C(idx, :);
noise = sqrt(sum((signal - centers).^2, 2));

noise_power = zeros(16, 1);
for k=1:16
    noise_power(k) = mean(noise(idx==k));
end

snr_point_cloud = sqrt(sum(C.^2, 2))./noise_power;
snr_total = mean(sqrt(sum(centers.^2, 2)))/mean(noise);

%% Save results for later use
save ssf_receive_signal_5.mat

