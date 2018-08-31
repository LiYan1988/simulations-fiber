function [data_mod_t_tmp, param] = generate_16qam(param, c)
% Generate 16QAM channel

% Generate bit stream, 0 at the start and end of the bit stream
data_bit = randi([0, 1], param.symbol_number(c)-2, param.bit_per_symbol(c));
%     data_bit = [zeros(1, param.bit_per_symbol(c)); data_bit; zeros(1, param.bit_per_symbol(c))];
param.data_bit_channel{c} = data_bit;

% Convert bits to integers and then to symbols
data_symbol = bi2de(data_bit);
param.data_symbol_channel{c} = data_symbol; % these are integers

data_mod_t_tmp = qammod(data_symbol, param.constellation_size(c));
data_mod_t_tmp = modnorm(data_mod_t_tmp, 'avpow', 1)*data_mod_t_tmp;

% The initial symbols in the constellation diagram
param.data_mod_symbol_channel{c} = data_mod_t_tmp;

% pad zeros at head and tail, prevent delay of filter distort signal
data_mod_t_tmp = [0; data_mod_t_tmp; 0];

% if it is 16QAM, use square-root RRC
param.filter_tx_channel{c} = rcosdesign(param.filter_parameter(c), ...
    param.symbol_in_filter(c), param.sample_per_symbol(c), 'sqrt');

% Delay of filter
delay_filter = (length(param.filter_tx_channel{c})-1)/2;
param.delay_filter_channel(c) = delay_filter;

% pass through filter
data_mod_t_tmp = upfirdn(data_mod_t_tmp, param.filter_tx_channel{c}, ...
    param.sample_per_symbol(c));

% remove head and tail
A = param.symbol_in_filter(c)/2; % half filter length
u = ((A-0.5)*param.sample_per_symbol(c)+1):...
    ((A-0.5)*param.sample_per_symbol(c)+param.fn); % index of the signal
data_mod_t_tmp = data_mod_t_tmp(u);

