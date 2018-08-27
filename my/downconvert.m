function [xt, xp, padding_length] = downconvert(signal, param, cidx)
% Downconvert cidx-th channel to baseband and low-pass filter
% Return baseband signal, xt, and samples xp

xt = signal.*exp(1i*2*pi*param.center_frequency_channel(cidx).*param.t);

% Pass through matching filter
xt = upfirdn(xt, param.filter_tx_channel{cidx}, 1, 1);
% Remove head and tail of the signal
xt = xt(param.filter_delay(cidx)+1:end-param.filter_delay(cidx));

% Zero padding to multiples of sampling rate
padding_length = ceil(length(xt)/param.sample_per_symbol(cidx))*...
    param.sample_per_symbol(cidx)-length(xt);
% reshape, each row is samples at one time shift
xp = reshape([xt; zeros(padding_length, 1)], ...
    param.sample_per_symbol(cidx), []);
