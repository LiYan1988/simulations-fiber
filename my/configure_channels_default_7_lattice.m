function param = configure_channels_default_7_lattice(param, ...
    power_dbm_qam, power_dbm_ook, ...
    bw_hz_qam, bw_hz_ook, grid_hz, N)
% Change channel parameters
%   1. power of 16QAM and OOK
%       power_dbm_16qam
%       power_dbm_ook
%   2. bandwidth of QAM and OOK
%       bw_ghz_ook
%       bw_ghz_qam
%   3. channel spacing
%   4. lattice channels
%       OOK-QAM-OOK-QAM-...-OOK-QAM-OOK
%       N is an odd number

% [Hz], spectrum grid size
spectrum_grid_size = grid_hz;

% channel type
channel_type = [repmat({'ook'}, 1, (N-1)/2); ...
    repmat({'16qam'}, 1, (N-1)/2)];
channel_type = channel_type(:);
channel_type(end+1) = {'ook'};

% [W], power of channel in time domain, in contrast to the frequency domain
% PSD measured in W/Hz
power_dbm = power_dbm_ook*ones(N, 1); % power of OOK
power_dbm(2:2:end) = power_dbm_qam;

% filter parameter
filter_parameter = 0.5*ones(1, N);
% For 16QAM use square-root RRC, then specify the roll-off factor
filter_parameter(2:2:end) = 0.2;

% symbol in filter
symbol_in_filter = 10*ones(1, N);

param = configure_channels_2(param, N, spectrum_grid_size, ...
    channel_type, power_dbm, filter_parameter, symbol_in_filter, ...
    bw_hz_ook, bw_hz_qam);