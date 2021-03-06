clc;
clear;
close all;

% OOK-QAM-OOK, power of OOK = -5 dBm

%% Fiber Parameters
% -------------- Primary parameters
param.fmax = 2*pi*1e9*(320); % [Hz], should be common multiples of all channels' bandwidths
param.fn = 2^17; % number of spectrum points

param.span_length = 82; % [km], span length
param.beta2 = -2.1683e-23; % [s^2/km], GVD, D=17 [ps/ns/km]
% S = 0.06*1e6; % [s/(m^2*km)], third order dispersion
param.wavelength = 1550*1e-9; % [m], reference wavelength
param.light_speed = 2.99792458*1e8; % [m/s], speed of light in vaccum
S = 4*pi*param.light_speed/param.wavelength^3*param.beta2;
param.beta3 = (S-4*pi*param.light_speed/param.wavelength^3*param.beta2)*...
    (param.wavelength^2/(2*pi*param.light_speed))^2; % [s^3/km], third order dispersion
param.gamma = 1.4; % [(W*km)^-1], nonlinear coefficient of SMF
param.alpha = log(10)*0.2/10; % [1/km] in linear, 0.2 dB/km, positive number
param.dz = 0.82; % [km]
param.ase_exist = true;
param.nsp = 1.2559; % [1] spontaneous emission factor, NF=5.563, nsp = 10^(NF/10)/2
param.h = 6.626*1e-34; % [J*s], [W*Hz^-2] Plank's constant
param.nu = param.light_speed/param.wavelength; % [Hz], light speed is in fiber, so 1.5 can bring it back to normal light speed

param.fbg_length_1 = 82; % [km], length of the first type DCF
param.fbg_length_2 = 82; % [km], length of the second type DCF

param.random_seed = 1;

%% Test
power_dbm_ook = -20:2:10;
bw_hz_ook = 10*1e9;
bw_hz_qam = 30*1e9;
grid_hz = 50*1e9;

time_elapsed = zeros(size(power_dbm_ook));
param_mp = cell(size(power_dbm_ook));


for n=1:length(power_dbm_ook)
    fprintf('%d of %d iterations\n', n, length(power_dbm_ook))
%     t = tic;
    param_tmp = configure_channels_default_8(param, ...
        power_dbm_ook(n), ...
        bw_hz_qam, ...
        bw_hz_ook, ...
        grid_hz, 1);
    
    % Generate Signal
    param_tmp = generate_signals(param_tmp);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
%     % plot sepctrum
%     figure; hold on; box on; grid on;
%     plot(param_tmp.f_plot, 10*log10(abs(param_tmp.data_mod_f_in).^2*1e12))
%     fprintf('sample per symbol %d, symbol number %d\n', ...
%         param_tmp.sample_per_symbol, param_tmp.symbol_number)
%     
%     % eye diagram
%     figure; hold on; box on; grid on;
%     plot(rem(param_tmp.t, param_tmp.sample_per_symbol*param_tmp.dt), ...
%         abs(param_tmp.data_mod_t_in).^2, '.')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Propagation through a link
    param_tmp = simulate_link2(param_tmp, 4);
    param_mp{n} = param_tmp;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % plot output spectrum
%     figure; hold on; box on; grid on;
%     plot(param_tmp.f_plot, 10*log10(abs(param_tmp.data_mod_f_current).^2*1e12))
%     fprintf('sample per symbol %d, symbol number %d\n', ...
%         param_tmp.sample_per_symbol, param_tmp.symbol_number)
%     
%     % received constellation
%     figure;
%     box on;
%     grid on;
%     hold on;
%     plot(param_tmp.tx_best_match{1}, '.')
%     plot(param_tmp.signal_tx_complex{1}+1e-9*1i, 'rx', 'linewidth', 2)
%     lim = max(max(abs(get(gca, 'xlim'))), max(abs(get(gca, 'ylim'))));
%     ylim([-lim, lim])
%     xlim([-lim, lim])
%     pbaspect([1 1 1])
%     
%     
%     % plot output eye diagram
%     figure; hold on; box on; grid on;
%     plot(rem(param_tmp.t, param_tmp.sample_per_symbol*param_tmp.dt), ...
%         abs(param_tmp.data_mod_t_dc).^2, '.')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

filename = 'ssf_a_2_compare_with_diego_ook_3_07.mat';
save(filename, '-v7.3')

%%
load(filename)
evm = zeros(size(power_dbm_ook));
snr = zeros(size(power_dbm_ook));
q1 = zeros(size(power_dbm_ook));
q2 = zeros(size(power_dbm_ook));
for n=1:length(power_dbm_ook)
    evm(n) = param_mp{n}.evm_channel;
    snr(n) = param_mp{n}.snr_channel;
    q1(n) = param_mp{n}.q_channel_1;
    q2(n) = param_mp{n}.q_channel_2;
end

figure;
hold on;
grid on;
box on;
% plot(power_dbm_ook, 10*log10(evm), 'o', 'linewidth', 2, 'displayname', 'EVM')
plot(power_dbm_ook, 10*log10(snr), 'o', 'linewidth', 2, 'displayname', 'SNR')
xlabel('Power (dBm)')
ylabel('SNR (dB)')
legend()

% figure;
% hold on;
% grid on;
% box on;
% plot(power_dbm_ook, q1, 'o', 'linewidth', 2, 'displayname', 'Q1')
% plot(power_dbm_ook, q2, 'o', 'linewidth', 2, 'displayname', 'Q2')
% legend()

% figure;
% hold on;
% box on;
% grid on;
% plot(q, 10*log10(evm), 'o', 'linewidth', 2, 'displayname', 'EVM vs Q')
% plot(q, snr, 'o', 'linewidth', 2, 'displayname', 'SNR vs Q')
% xlabel('Q')
% ylabel('EVM (dB)')
% l = legend();
% l.Location = 'northwest';
%% 
figure;
hold on;
box on;
grid on;
plot(param_mp{1}.f_plot, 10*log10(abs(param_mp{1}.data_mod_f_in).^2*1e12))
% plot(param_mp{1}.f_plot, smooth(10*log10(abs(param_mp{1}.data_mod_f_in).^2*1e12), 1000), 'linewidth', 2)
% set(gca, 'yscale', 'log')
grid on;
