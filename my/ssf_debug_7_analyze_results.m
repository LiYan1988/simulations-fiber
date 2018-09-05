clc;
clear;
close all;

% Analyze results of uniform power

%% 
load variable_power_11channels.mat

%%
snr_16qam = zeros(1, length(power_dbm));
snr_ook = zeros(1, length(power_dbm));
for n=1:length(power_dbm)
    snr_16qam(n) = param_mp{n}.snr_channel(6);
    snr_ook(n) = param_mp{n}.snr_channel(5);
%     x = param_mp{n}.signal_received_constellation_derotate{6};
%     figure;
%     plot(x(:, 1), x(:, 2), '.')
end

figure;
hold on;
box on;
grid on;
plot(power_dbm, 10*log10(snr_16qam), 'displayname', '16QAM')
plot(power_dbm, 10*log10(snr_ook), 'displayname', 'OOK')
legend('show')
xlabel('Power (dBm)')
ylabel('SNR (dB)')
%%
param = param_mp{5}; % this is the optimal power

% plot spectrum
figure;
plot(param.f_plot, 10*log10(abs(param.data_mod_f_current).^2))

% constellation of OOK 
ook = param.signal_received_constellation_derotate{5};
ook_centers = param.cloud_centers_derotation{5};
figure;
hold on;
plot(ook(:, 1), ook(:, 2), '.')
plot(ook_centers(:, 1), ook_centers(:, 2), 'x', 'linewidth', 2)

% constellation of 16QAM
qam = param.signal_received_constellation_derotate{6};
qam_centers = param.cloud_centers_derotation{6};
figure;
hold on;
plot(qam(:, 1), qam(:, 2), '.')
plot(qam_centers(:, 1), qam_centers(:, 2), 'x', 'linewidth', 2)

%% 
