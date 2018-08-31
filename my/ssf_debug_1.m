clc;
clear;
close all;

% Simulate varying channel spacing, all channel power equals -1 dBm

%% Fiber Parameters
% -------------- Primary parameters
param.fmax = 2*pi*640*1e9; % [Hz], should be common multiples of all channels' bandwidths
param.fn = 2^17; % number of spectrum points

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
param.ase_exist = true;
param.nsp = 1.8; % [1] spontaneous emission factor, NF=5.5
param.h = 6.626*1e-34; % [J*s], [W*Hz^-2] Plank's constant
param.nu = param.light_speed/param.wavelength*1.5; % [Hz], light speed is in fiber, so 1.5 can bring it back to normal light speed

%% Channel Parameters
% Channel specific parameters, n channels should have n sets of parameters
N = 3; % number of channels, should be an odd number
param.bandwidth_channel = 10*1e9*ones(1, N); % [Hz]
param.bandwidth_channel((N-1)/2+1) = 32*1e9; % CUT


%%%%%%%%% constellation_size is not good, will not use it %%%%%%%%%%
param.constellation_size = 2*ones(1, N); % 2=BPSK; 4=QPSK; 8=8QAM; 16=16QAM; 32=32QAM; 64=64QAM
% param.constellation_size = [2, 4, 8, 16, 32];
param.constellation_size((N-1)/2+1) = 16;
%%%%%%%%% should directlly specify modulation format by string %%%%%%%%%%
param.channel_type = [repmat({'ook'}, (N-1)/2, 1); {'16qam'}; ...
    repmat({'ook'}, (N-1)/2, 1)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



param.spectrum_grid_size = 50*1e9; % [Hz], spectrum grid size
param.center_frequency_channel = param.spectrum_grid_size*(linspace(0, N-1, N)-(N-1)/2);
param.power_channel_time = 10^(-1/10)/1e3*ones(N, 1); % [W], power of channel in time domain, in contrast to the frequency domain PSD measured in W/Hz

%%%%%%%%%%%%%%%%%%%% Change the filter parameter configuration %%%%%%%%%
% For OOK channels, specify the bandwidth-symbol time product
param.filter_parameter = 1*ones(1, N);
% For 16QAM use square-root RRC, then specify the roll-off factor
param.filter_parameter((N-1)/2+1) = 0.2; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


param.symbol_in_filter = 20*ones(1, N); % length of impulse response in symbol

%%%%%%%%%%%%%% square-root RRC is always used for 16QAM, no need to specify
% param.shape_filter = cell(1, N); % shape of raised
% param.shape_filter(:) = {'sqrt'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Random seed
param.random_seed = 2394759; % input to rng

%% Vary channel spacing
spacing_step = ((-10:2:30)+50)*1e9;
param_mp = cell(1, length(spacing_step)); % [dBm], power of each channel

param_temp = param;

% Generate Signal
param_temp = generate_signals(param_temp);

% Propagation through a link
param_temp = simulate_link1(param_temp);

param_mp{k} = param_temp;

%% Save results
save mp11_simulation_varying_channel_spacing.mat

%% Plot results
% load mp8_simulation_varying_16qam_power.mat
n_mp = length(param_mp);
cidx = (param_mp{1}.channel_number+1)/2;
snr_16qam = zeros(n_mp, 1);
snr_5ook = zeros(n_mp, 1);
for n=1:n_mp
    snr_16qam(n) = param_mp{n}.snr_total{cidx}(1);
    snr_5ook(n) = param_mp{n}.snr_total{cidx-1}(1);
end

figure;
box on;
grid on;
hold on;
title('OOK Power @-1dBm')
h1 = plot(spacing_step/1e9, 10*log10(snr_16qam), 'displayname', '16 QAM', ...
    'linewidth', 2);
h2 = plot(spacing_step/1e9, 10*log10(snr_5ook), 'displayname', '5th OOK', ...
    'linewidth', 2);
xlabel('Channel spacing (GHz)')
ylabel('SNR (dB)')
legend([h1, h2])
pbaspect([7 4 1])