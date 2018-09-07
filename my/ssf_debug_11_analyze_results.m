clc;
clear;
close all;

% Simulate varying channel number, including even number of channels
% 16QAM power 0dBm, OOK power -6dBm

%% Load results
load debug_11_variable_channel_number.mat

%%
snr_qam = zeros(size(param_mp));
ser_qam = zeros(size(param_mp));
snr_ook = zeros(size(param_mp));

for n=1:length(param_mp)
    k = (n+2-mod(n, 2))/2;
    snr_qam(n) = param_mp{n}.snr_channel(k);
    ser_qam(n) = param_mp{n}.ser_channel(k);
    if n>1
        snr_ook(n) = param_mp{n}.snr_channel(k-1);
        ser_ook(n) = param_mp{n}.ser_channel(k-1);
    end
end

figure; hold on; box on; grid on;
plot(snr_qam, 'displayname', '16QAM', 'linewidth', 2)
plot(snr_ook, 'displayname', 'OOK', 'linewidth', 2)
xlabel('Number of channels')
ylabel('Linear SNR')
legend()
xlim([1, 40])
% figure; hold on; box on; grid on;
% plot(ser_qam)
% plot(ser_ook)

figure; hold on; box on; grid on;
plot(1./snr_qam, 'displayname', '16QAM', 'linewidth', 2)
plot(1./snr_ook, 'displayname', 'OOK', 'linewidth', 2)
xlabel('Number of channels')
ylabel('Inverse SNR')
legend('Location', 'northwest')
xlim([1, 40])
title('1/SNR')