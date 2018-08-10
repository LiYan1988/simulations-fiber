function [xt_dc, xf_dc, xt] = dispersion_compensation(signal, param, cidx)
% Compensate dispersion for cidx-th channel
% xt_dc: time domain signal after compensation
% xf_dc: frequency domain signal after compensation
% xt: time domain signal before compensation
[xt, ~, ~] = downconvert(signal, param, cidx);

% a all-pass filter compensating for dispersion
dc_filter = exp(-0.5*1i*param.beta2*(param.f+...
    2*pi*param.center_frequency_channel(cidx)).^2*param.span_length); 
xf_dc = ft(xt, param.df);
xf_dc = xf_dc.*dc_filter;
xt_dc = ift(xf_dc, param.df);