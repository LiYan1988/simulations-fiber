clc;
clear;
close all;

% OOK-QAM-OOK, power of OOK = -5 dBm

%% Fiber Parameters
% -------------- Primary parameters
param.fmax = 2*pi*1e9*119; % [Hz], should be common multiples of all channels' bandwidths
param.fn = 2^15; % number of spectrum points

param.fmax = 2*pi*1e9*(60-1); % [Hz], should be common multiples of all channels' bandwidths
param.fn = 2^14; % number of spectrum points

param.span_length = 82; % [km], span length
param.beta2 = -2.1683e-23; % [s^2/km], GVD, D=17 [ps/ns/km]
% S = 0.06*1e6; % [s/(m^2*km)], third order dispersion
param.wavelength = 1550*1e-9; % [m], reference wavelength
param.light_speed = 2.99792458*1e8; % [m/s], speed of light in fiber
S = 4*pi*param.light_speed/param.wavelength^3*param.beta2;
param.beta3 = (S-4*pi*param.light_speed/param.wavelength^3*param.beta2)*...
    (param.wavelength^2/(2*pi*param.light_speed))^2; % [s^3/km], third order dispersion
% d3=1i*beta3/6*(2*pi*FF).^3;
param.gamma = 1.4; % [(W*km)^-1], nonlinear coefficient of SMF
param.alpha = log(10)*0.2/10; % [1/km] in linear, 0.2 dB/km, positive number
param.dz = 0.82; % [km]
param.ase_exist = true;
param.nsp = 1.2559; % [1] spontaneous emission factor, NF=5.563, nsp = 10^(NF/10)/2
param.h = 6.626*1e-34; % [J*s], [W*Hz^-2] Plank's constant
param.nu = param.light_speed/param.wavelength; % [Hz], light speed is in fiber, so 1.5 can bring it back to normal light speed

param.fbg_length_1 = 82; % [km], length of the first type DCF
param.fbg_length_2 = 82; % [km], length of the second type DCF

param.random_seed = 1; %54790;

%% Test
power_dbm_ook = -5;
power_dbm_qam = -10:2:10;
bw_hz_ook = 10*1e9;
bw_hz_qam = 30*1e9;
grid_hz = 25*1e9;

time_elapsed = zeros(size(power_dbm_qam));
param_mp = cell(size(power_dbm_qam));

parfor n=1:length(power_dbm_qam)
    fprintf('%d of %d iterations\n', n, length(power_dbm_qam))
%     t = tic;
    param_tmp = configure_channels_default_6(param, ...
        power_dbm_qam(n), ...
        power_dbm_ook, ...
        bw_hz_qam, ...
        bw_hz_ook, ...
        grid_hz, 1);
    
    % Generate Signal
    param_tmp = generate_signals(param_tmp);
    
    % plot sepctrum
%     figure; hold on; box on; grid on;
%     plot(param_tmp.f_plot, 10*log10(abs(param_tmp.data_mod_f_in).^2*1e12))
%     set(gca, 'yscale', 'log')
    
    % Propagation through a link
    param_tmp = simulate_link2(param_tmp, 4);
%     time_elapsed(n) = toc(t);
%     fprintf('Running time %.f minutes, total time %.f minutes\n', ...
%         time_elapsed(n)/60, sum(time_elapsed)/60)
%     fprintf('-------------------------------------------\n')
    param_mp{n} = param_tmp;
    
    % plot
%     plot(param.f_plot, abs(param.data_mod_f_in).^2)
%     set(gca, 'yscale', 'log')
%     
%     scatterplot(param.signal_rx_complex{1})
%     evm(n) = param.evm_channel;
end

%%
evm = zeros(size(power_dbm_qam));
snr = zeros(size(power_dbm_qam));
for n=1:length(power_dbm_qam)
    evm(n) = param_mp{n}.evm_channel;
    snr(n) = param_mp{n}.snr_channel;
end

figure;
hold on;
grid on;
box on;
%plot(power_dbm_qam, 10*log10(evm), 'o', 'linewidth', 2, 'displayname', 'EVM')
plot(power_dbm_qam, 10*log10(snr), 'o', 'linewidth', 2, 'displayname', 'SNR')
xlabel('Power (dBm)')
ylabel('EVM (dB)')
legend()

%% 
figure;
box on;
grid on;
plot(param_mp{1}.f_plot, abs(param_mp{1}.data_mod_f_in).^2)
set(gca, 'yscale', 'log')
grid on;
