clc;
clear;
close all;

% Prepare to do more simulations with multiple parameters
% Parameters to be changed:
% 1. number of channels
% 2. channel spacing
% 3. channel power
%   a. uniform power, but changing
%   b. especially power of OOK to coherent channels, i.e., different power
%   between OOk and 16QAM
% 4. N of spans
% Don't worry about dispersion or nonlinear coefficient, those are
% relatively fixed.

% In this file, change relative power of OOK to coherent channels
% The power of the 16QAM channel is fixed at the optimal power (-1 dBm)

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
param.power_channel_time = 10^(-3/10)/1e3*ones(N, 1); % [W], power of channel in time domain, in contrast to the frequency domain PSD measured in W/Hz

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
    % Fix the power of 16QAM to -3 dBm
    param_temp.power_channel_time((N-1)/2+1) = 10^(-1/10)/1e3*ones(N, 1); % [W]
    
    % Generate Signal
    param_temp = generate_signals(param_temp);
    
    % Propagation through a link
    param_temp = simulate_link1(param_temp);
    
    param_mp{k} = param_temp;
end

%% Save results
save simulation_uniform_power_1.mat

%% Plot results
n_mp = length(param_mp);
cidx = (param_mp{1}.channel_number+1)/2;
snr_total = zeros(n_mp, 1);
for n=1:n_mp
    snr_total(n) = param_mp{n}.snr_total{cidx}(1);
end

figure;
title('SNR of 16QAM (all channels have the same power)')
plot(power_step, 10*log10(snr_total))
xlabel('Power (dBm)')
ylabel('SNR with only NLI (dB)')