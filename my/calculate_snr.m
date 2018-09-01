function param = calculate_snr(param)
% calculate snr for each channel

% channel SNRs
param.snr_channel = zeros(1, param.channel_number);
% channel linear signal powers
param.power_signal_linear = zeros(1, param.channel_number);
% channel linear noise power
param.power_noise_linear = zeros(1, param.channel_number);

for cidx = 1:param.channel_number
    % received symbols
    rx = param.signal_rx_complex{cidx};
    % transmitted symbols
    tx = param.signal_tx_complex{cidx};
    % unique transmitted symbols
    tx_unique = unique(tx);
    
    power_signal = 0; % signal power
    power_noise = 0; % noise power
    % iterate over all unique transmitted symbols
    for n = 1:length(tx_unique)
        tmp_signal = rx(tx==tx_unique(n));
        tmp_ratio = length(tmp_signal)/length(rx);
        power_signal = power_signal + abs(mean(tmp_signal))^2*tmp_ratio;
        power_noise = power_noise + abs(std(tmp_signal))^2*tmp_ratio;
    end
    param.snr_channel(cidx) = power_signal/power_noise;
    param.power_signal_linear(cidx) = power_signal;
    param.power_noise_linear(cidx) = power_noise;
end