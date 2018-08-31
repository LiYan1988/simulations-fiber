function [xt] = downconvert(signal, param, cidx)
% Downconvert cidx-th channel to baseband and low-pass filter
% Return baseband signal, xt, and samples xp

xt = signal.*exp(1i*2*pi*param.center_frequency_channel(cidx).*param.t);

% Pass through matching filter
xt = upfirdn(xt, param.filter_tx_channel{cidx}, 1, 1);
% Remove head and tail of the signal
s = (param.filter_delay(cidx)+1):(param.filter_delay(cidx)+param.fn);
xt = xt(param.filter_delay(cidx)+1:end-param.filter_delay(cidx));
