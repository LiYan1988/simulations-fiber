clc;
clear;
close all;

% analyze variable channel number

%% Load results
load 'debug_10_variable_channel_number.mat'

%%
number_of_channel = 1:2:21; % total number of channels
snr_qam = zeros(size(number_of_channel));
snr_ook = zeros(size(number_of_channel));

for n=1:length(number_of_channel)
    snr_qam(n) = param_mp{n}.snr_channel(n);
    if number_of_channel(n)>1
        snr_ook(n) = param_mp{n}.snr_channel(n-1);
    end
end

figure;
hold on;
box on;
grid on;
plot(number_of_channel, snr_qam, 'displayname', '16QAM')
plot(number_of_channel, snr_ook, 'displayname', 'OOK')
xlabel('Number of channels')
ylabel('Linear SNR')
legend()
xlim([1, 21])
xticks(number_of_channel)