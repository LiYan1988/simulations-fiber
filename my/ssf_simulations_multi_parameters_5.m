clc;
clear;
close all;

% Simulate uniform power

%% Fiber Parameters
% -------------- Primary parameters
param.fmax = 2*pi*500*1e9; % [Hz]
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
N = 11; % number of channels, should be an odd number
param.bandwidth_channel = 10*1e9*ones(1, N); % [Hz]
param.bandwidth_channel((N-1)/2+1) = 32*1e9; % CUT
param.constellation_size = 2*ones(1, N); % 2=BPSK; 4=QPSK; 8=8QAM; 16=16QAM; 32=32QAM; 64=64QAM
% param.constellation_size = [2, 4, 8, 16, 32];
param.constellation_size((N-1)/2+1) = 16;
param.spectrum_grid_size = 50*1e9; % [Hz], spectrum grid size
param.center_frequency_channel = param.spectrum_grid_size*(linspace(0, N-1, N)-(N-1)/2);
param.power_channel_time = 10^(-1/10)/1e3*ones(N, 1); % [W], power of channel in time domain, in contrast to the frequency domain PSD measured in W/Hz

% Filter parameters of each channel, assume raised cosine filters
param.roll_off_filter = 0.2*ones(1, N); % Roll-off factor of raised cosine filter
param.symbol_in_filter = 100*ones(1, N); % length of impulse response in symbol
param.shape_filter = cell(1, N); % shape of raised
param.shape_filter(:) = {'sqrt'};

% Random seed
param.random_seed = 2394759; % input to rng

%% Vary uniform power of all the channels
power_step = -10:1:10;
param_mp = cell(1, length(power_step)); % [dBm], power of each channel

parfor k=1:length(power_step)
    param_temp = param;
    
    % Change the uniform power
    param_temp.power_channel_time = 10^(power_step(k)/10)/1e3*ones(N, 1); % [W]
%     % Fix the power of 16QAM to -1 dBm
%     param_temp.power_channel_time((N-1)/2+1) = 10^(-1/10)/1e3; % [W]
    
    % Generate Signal
    param_temp = generate_signals(param_temp);
    
    % Propagation through a link
    param_temp = simulate_link1(param_temp);
    
    param_mp{k} = param_temp;
end

%% Save results
save simulation_uniform_power.mat

%% Plot results
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
title('Uniform Launch Power')
h1 = plot(power_step, 10*log10(snr_16qam), 'displayname', '16 QAM', ...
    'linewidth', 2);
h2 = plot(power_step, 10*log10(snr_5ook), 'displayname', '5th OOK', ...
    'linewidth', 2);
xlabel('Power (dBm)')
ylabel('SNR (dB)')
legend([h1, h2])
pbaspect([7 4 1])

%% Constellation at -1dBm
% 16QAM
constellation_16qam = param_mp{10}.signal_received_constellation_derotate{cidx};
constellation_16qam_centers = param_mp{10}.cloud_centers_derotation{cidx};
figure;
hold on;
box on;
grid on;
plot(constellation_16qam(:, 1), constellation_16qam(:, 2), '.')
plot(constellation_16qam_centers(:, 1), constellation_16qam_centers(:, 2), 'x', 'linewidth', 2)
xlabel('In-Phase')
ylabel('Quadrature')
xlim([-0.4, 0.4])
ylim([-0.4, 0.4])
pbaspect([1, 1, 1])


% 5th OOK
constellation_5ook = param_mp{10}.signal_received_constellation_derotate{cidx-1};
constellation_5ook_centers = param_mp{10}.cloud_centers_derotation{cidx-1};
figure;
hold on;
box on;
grid on;
plot(constellation_5ook(:, 1), constellation_5ook(:, 2), '.')
plot(constellation_5ook_centers(:, 1), constellation_5ook_centers(:, 2), 'x', 'linewidth', 2)
xlabel('In-Phase')
ylabel('Quadrature')
xlim([-0.4, 0.4])
ylim([-0.4, 0.4])
pbaspect([1, 1, 1])

%% Constellation with histogram color map
% 16 QAM
% [hist_n, hist_c] = hist3(constellation_16qam, 'Nbins', [100, 100]);
[hist_n, hist_c] = hist3(constellation_16qam, ...
    'Ctrs', {linspace(-0.4, 0.4, 100), linspace(-0.4, 0.4, 100)});
hist_n = hist_n/sum(hist_n(:));
% contour(hist_n)

imagesc(hist_c{1}, fliplr(hist_c{2}), hist_n)
colormap(hot)
xlabel('In-Phase')
ylabel('Quadrature')
pbaspect([1, 1, 1])

%%
% 5th OOK
[hist_n, hist_c] = hist3(constellation_5ook, ...
    'Ctrs', {linspace(-0.4, 0.4, 100), linspace(-0.4, 0.4, 100)});
hist_n = hist_n/sum(hist_n(:));
% contour(hist_n)

h = imagesc(hist_c{1}, -hist_c{2}, hist_n');
colormap(hot)
xlabel('In-Phase')
ylabel('Quadrature')
pbaspect([1, 1, 1])