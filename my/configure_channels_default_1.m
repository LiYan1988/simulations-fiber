function param = configure_channels_default_1(param, N)
% Only change the number of channels, other parameters are the default
% values
% The optimal power is used
% OOK power is -6 dBm
% 16QAM power is 0 dBm

% [Hz], spectrum grid size
spectrum_grid_size = 50*1e9;

% channel type
if mod(N, 2)==1
    channel_type = [repmat({'ook'}, (N-1)/2, 1); {'16qam'}; ...
        repmat({'ook'}, (N-1)/2, 1)];
elseif mod(N, 2)==0
    channel_type = [repmat({'ook'}, N/2, 1); {'16qam'}; ...
        repmat({'ook'}, N/2-1, 1)];
end

% [W], power of channel in time domain, in contrast to the frequency domain
% PSD measured in W/Hz
power_dbm = -6*ones(N, 1);
if mod(N, 2)==1
    power_dbm((N-1)/2+1) = 0;
elseif mod(N, 2)==0
    power_dbm(N/2+1) = 0;
end

% filter parameter
filter_parameter = 0.7*ones(1, N);
% For 16QAM use square-root RRC, then specify the roll-off factor
if mod(N, 2)==1
    filter_parameter((N-1)/2+1) = 0.2;
elseif mod(N, 2)==0
    filter_parameter(N/2+1) = 0.2;
end

% symbol in filter
symbol_in_filter = 10*ones(1, N);

param = configure_channels(param, N, spectrum_grid_size, ...
    channel_type, power_dbm, filter_parameter, symbol_in_filter);