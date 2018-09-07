function param = configure_channels_default_5(param, ...
    power_dbm_qam, power_dbm_ook, ...
    bw_ghz_qam, bw_ghz_ook, grid_ghz)
% Change channel parameters
%   1. power of 16QAM and OOK
%       power_dbm_16qam
%       power_dbm_ook
%   2. bandwidth of QAM and OOK
%       bw_ghz_ook
%       bw_ghz_qam
%   3. channel spacing

N = 3;

% [Hz], spectrum grid size
spectrum_grid_size = grid_ghz;

% channel type
channel_type = [repmat({'ook'}, (N-1)/2, 1); {'16qam'}; ...
    repmat({'ook'}, (N-1)/2, 1)];

% [W], power of channel in time domain, in contrast to the frequency domain
% PSD measured in W/Hz
power_dbm = power_dbm_ook*ones(N, 1); % power of OOK
power_dbm((N-1)/2+1) = power_dbm_qam;

% filter parameter
filter_parameter = 0.7*ones(1, N);
% For 16QAM use square-root RRC, then specify the roll-off factor
filter_parameter((N-1)/2+1) = 0.2;

% symbol in filter
symbol_in_filter = 10*ones(1, N);

param = configure_channels_2(param, N, spectrum_grid_size, ...
    channel_type, power_dbm, filter_parameter, symbol_in_filter, ...
    bw_ghz_ook, bw_ghz_qam);