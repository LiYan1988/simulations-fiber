function param = configure_channels_default_8(param, ...
    power_dbm_ook, bw_hz_qam, bw_hz_ook, grid_hz, N)
% Change channel parameters
%   1. power of 16QAM and OOK
%       power_dbm_16qam
%       power_dbm_ook
%   2. bandwidth of QAM and OOK
%       bw_ghz_ook
%       bw_ghz_qam
%   3. channel spacing
%   4. number of channels
%       N is an odd number
%   5. only OOK channels

% [Hz], spectrum grid size
spectrum_grid_size = grid_hz;

% channel type
channel_type = [repmat({'ook'}, (N-1)/2, 1); {'ook'}; ...
    repmat({'ook'}, (N-1)/2, 1)];

% [W], power of channel in time domain, in contrast to the frequency domain
% PSD measured in W/Hz
power_dbm = power_dbm_ook*ones(N, 1); % power of OOK

% filter parameter
filter_parameter = 0.5*ones(1, N);

% symbol in filter
symbol_in_filter = 10*ones(1, N);

param = configure_channels_2(param, N, spectrum_grid_size, ...
    channel_type, power_dbm, filter_parameter, symbol_in_filter, ...
    bw_hz_ook, bw_hz_qam);