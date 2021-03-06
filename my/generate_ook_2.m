function [data_mod_t_tmp, param] = generate_ook_2(param, c)
% Generate NRZ OOK signal

% Generate bit stream, 0 at the start and end of the bit stream
data_bit = randi([0, 1], param.symbol_number(c)-2, param.bit_per_symbol(c));
data_bit_tmp = [0; data_bit; 0];
% data_mod_t_tmp = ones(param.symbol_number(c), param.sample_per_symbol(c));
% mask = repmat(data_bit_tmp, 1, param.sample_per_symbol(c));
% data_mod_t_tmp = (data_mod_t_tmp.*mask)';
% data_mod_t_tmp = data_mod_t_tmp(:);

% NRZ OOK without carrier suppress
data_mod_t_tmp = repmat(data_bit_tmp, 1, param.sample_per_symbol(c))';
% data_mod_t_tmp(2:end, :) = 0;
% data_mod_t_tmp = data_mod_t_tmp(:);


% if it is OOK, use Gauss filter 
param.filter_tx_channel{c} = gaussdesign(param.filter_parameter(c), ...
    param.symbol_in_filter(c), param.sample_per_symbol(c));

% Delay of filter
delay_filter = (length(param.filter_tx_channel{c})-1)/2;
param.delay_filter_channel(c) = delay_filter;

% pass through low-pass filter
data_mod_t_tmp = upfirdn(data_mod_t_tmp, param.filter_tx_channel{c});

% remove head and tail
s = (param.delay_filter_channel(c)+1):...
    (param.delay_filter_channel(c)+length(param.t));
data_mod_t_tmp = data_mod_t_tmp(s);

data_mod_t_tmp = 0.5*(cos(pi*data_mod_t_tmp)).';

% suppress carrier
% xf = ft(data_mod_t_tmp, param.df);
% xf(length(xf)/2+1) = xf(length(xf)/2+1)*0.01;
% data_mod_t_tmp = ift(xf, param.df);

% record generated bits and symbols
param.data_bit_channel{c} = data_bit;
param.data_mod_symbol_channel{c} = data_bit; % these are constellation points
param.data_symbol_channel{c} = data_bit; % these are integers