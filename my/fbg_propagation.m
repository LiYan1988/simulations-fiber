function [signal_t_out, signal_f_out] = fbg_propagation(...
    signal, beta2, beta3, span_length, param)
% pass signal through FBG
% input: beta2, beta3 are dispersion coefficient of the equivalent DCF,
% span_length is the equivalent length of DCF
% output: time and frequency domain signals

fbg_dispersion = exp(0.5*1i*beta2*param.f.^2*span_length).*...
    exp(1i/6*beta3*param.f.^3*span_length);

signal_t_out = ift(ft(signal, param.df).*fbg_dispersion, param.df);
signal_f_out = ft(signal_t_out, param.df);