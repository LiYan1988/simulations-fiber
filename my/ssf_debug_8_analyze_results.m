clc;
clear;
close all;

%% Load results
load variable_power_11channels_2.mat

%% Plot results
snr_qam = zeros(size(param_mp));
ser_qam = zeros(size(param_mp));
snr_ook = zeros(size(param_mp));
ser_ook = zeros(size(param_mp));

tic
for k=1:size(snr_qam, 1)
    fprintf('%d\n', k)
    for m=1:size(snr_qam, 2)
        param_mp{k,m} = calculate_ser_evm(param_mp{k,m});
        N = param_mp{k,m}.channel_number;
        snr_qam(k,m) = param_mp{k,m}.snr_channel((N-1)/2+1);
        ser_qam(k,m) = param_mp{k,m}.ser_channel((N-1)/2+1);
        snr_ook(k,m) = param_mp{k,m}.snr_channel((N-1)/2);
        ser_ook(k,m) = param_mp{k,m}.ser_channel((N-1)/2);
    end
end
time_elapsed = toc;

%%
save('variable_power_11channels_2.mat','-v7.3')

%% 
time_elapsed = zeros(10, 1);
parfor n=1:10
    t = tic;
    surf(peaks(40));
    time_elapsed(n) = toc(t);
end