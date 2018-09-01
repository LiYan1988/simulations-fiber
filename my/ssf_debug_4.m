clc;
clear;
close all;

load matlab2.mat

%%
% cidx = 2;
% rx = param.signal_rx_complex{cidx};
% tx = param.signal_tx_complex{cidx};
% tx_unique = unique(tx);
% 
% power_signal = 0; % power of signal
% power_noise = 0; % power of noise
% power_signal2 = 0;
% for n=1:length(tx_unique)
%     tmp_signal = rx(tx==tx_unique(n));
%     tmp_ratio = length(tmp_signal)/length(rx);
%     power_signal = power_signal + abs(mean(tmp_signal))^2*tmp_ratio;
%     power_noise = power_noise + abs(std(tmp_signal))^2*tmp_ratio;
% end
% 
% snr = power_signal/power_noise;

param = calculate_snr(param);

save matlab3.mat