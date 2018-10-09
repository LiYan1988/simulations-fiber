function [xt_dc, xf_dc, xt] = dispersion_compensation(signal, cidx, ...
    beta2, beta3, span_length, param)
% Compensate all the dispersion for cidx-th channel in SSMF
% beta2, beta3 are dispersion coefficients, length is the equivalent fiber
% length.
% NOTE: beta2 and beta3 are dispersion coefficients of SSMF, not DCF
% xt_dc: time domain signal after compensation
% xf_dc: frequency domain signal after compensation
% xt: time domain signal before compensation

% Downconvert signal to baseband
[xt] = downconvert(signal, param, cidx);

% For OOK channel, just take the absolute value, whereas for QAM channel,
% compensate residual dispersion
if false %strcmp(param.channel_type{cidx}, 'ook')
    xt_dc = abs(xt).^2;
    xf_dc = ft(xt_dc, param.df);
else
    % a all-pass filter compensating for 2nd and 3rd order dispersion
    % frequency axis
    tmp_f = param.f+2*pi*param.center_frequency_channel(cidx);
    dc_filter = exp(-0.5*1i*beta2*(tmp_f).^2*span_length).*...
        exp(-1i/6*beta3*(tmp_f).^3*span_length);
    
    xf_dc = ft(xt, param.df);
    xf_dc = xf_dc.*dc_filter;
    xt_dc = ift(xf_dc, param.df);
end