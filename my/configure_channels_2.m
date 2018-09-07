function param = configure_channels_2(param, N, spectrum_grid_size, ...
    channel_type, power_dbm, filter_parameter, symbol_in_filter,...
    ook_bw, qam_bw)
% Configure channel parameters

% Constellation size and type
% OOK: 2, 16QAM: 16
assert(length(channel_type)==N, ...
    sprintf('Length of channel_type should be %d', N))
param.channel_type = cell(1, N); % strings
param.constellation_size = zeros(1, N); % int
param.bandwidth_channel = zeros(1, N);
for n=1:N
    param.channel_type{n} = channel_type{n};
    if strcmp(channel_type{n}, 'ook')
        param.constellation_size(n) = 2;
        param.bandwidth_channel(n) = ook_bw;
    elseif strcmp(channel_type{n}, '16qam')
        param.constellation_size(n) = 16;
        param.bandwidth_channel(n) = qam_bw;
    end
end

% [Hz], spectrum grid size
param.spectrum_grid_size = spectrum_grid_size;
if mod(N, 2)==1
    param.center_frequency_channel = param.spectrum_grid_size*...
        (linspace(0, N-1, N)-(N-1)/2);
elseif mod(N, 2)==0
    param.center_frequency_channel = param.spectrum_grid_size*...
        (linspace(0, N-1, N)-N/2);
end

% [W], power of channel in time domain, in contrast to the frequency domain
% PSD measured in W/Hz
assert(length(power_dbm)==N, ...
    sprintf('Length of power_dbm should be %d', N))
param.power_channel_time = 10.^(power_dbm/10)/1e3;


% filter parameter
% Gaussian FIR for OOK
% square-root RRC for 16QAM
assert(length(filter_parameter)==N, ...
    sprintf('Length of filter_parameter should be %d', N))
param.filter_parameter = filter_parameter;

% number of symbols in filter
assert(length(symbol_in_filter)==N, ...
    sprintf('Length of symbol_in_filter should be %d', N))
param.symbol_in_filter = symbol_in_filter;

