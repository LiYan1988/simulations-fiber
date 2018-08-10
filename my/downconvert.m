function [xt, xp, padding_length] = downconvert(signal, param, cidx)
% Downconvert cidx-th channel to baseband and low-pass filter
% Return baseband signal, xt, and samples xp

% Downconvert
% if strcmp(mode, 'in')
%     x = param.data_mod_t_in.*exp(1i*2*pi*param.center_frequency_channel(cidx).*param.t);
% elseif strcmp(mode, 'current')
%     x = param.data_mod_t_current.*exp(1i*2*pi*param.center_frequency_channel(cidx).*param.t);
% end
x = signal.*exp(1i*2*pi*param.center_frequency_channel(cidx).*param.t);

% Spectrum domain
xf = ft(x, param.df);
% Low-pass filter with an ideal rectangular filter
xf = xf.*param.f_mask;
% Back to time domain
xt = ift(xf, param.df);

% Zero padding to multiples of sampling rate
padding_length = ceil(length(xt)/param.sample_per_symbol(cidx))*...
    param.sample_per_symbol(cidx)-length(xt);
% reshape, each row is samples at one time shift
xp = reshape([xt; zeros(padding_length, 1)], ...
    param.sample_per_symbol(cidx), []);
