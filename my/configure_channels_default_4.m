function param = configure_channels_default_4(param, N)
% Only change the number of channels, other parameters are the default
% values
% N is an even number
% The optimal power is used
% OOK power is -6 dBm
% 16QAM power is 0 dBm

% [Hz], spectrum grid size
spectrum_grid_size = 50*1e9;

% channel type
channel_type = [repmat({'ook'}, N/2, 1); {'16qam'}; ...
    repmat({'ook'}, N/2-1, 1)];

% [W], power of channel in time domain, in contrast to the frequency domain
% PSD measured in W/Hz
power_dbm = -6*ones(N, 1); 
power_dbm(N/2+1) = 0;

% filter parameter
filter_parameter = 0.7*ones(1, N);
% For 16QAM use square-root RRC, then specify the roll-off factor
filter_parameter(N/2+1) = 0.2;

% symbol in filter
symbol_in_filter = 10*ones(1, N);

param = configure_channels(param, N, spectrum_grid_size, ...
    channel_type, power_dbm, filter_parameter, symbol_in_filter);