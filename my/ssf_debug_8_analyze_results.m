clc;
clear;
close all;

%% Load results
load variable_power_11channels_2.mat

%% Plot results
snr_qam = zeros(size(param_mp));
snr_ook = zeros(size(param_mp));

for k=1:size(snr_qam, 1)
    for m=1:size(snr_qam, 2)
        N = param_mp{k,m}.channel_number;
        snr_qam(k,m) = param_mp{k,m}.snr_channel((N-1)/2+1);
        snr_ook(k,m) = param_mp{k,m}.snr_channel((N-1)/2);
    end
end