function fwhm = calculate_fir_bandwidth(param, cidx)
% calculate FIR Full width at half maximum (FWHM) bandwidth
% fwhm is in Hz
fir = param.filter_tx_channel{cidx};
[h, w] = freqz(fir);
h = abs(h);
h = h/max(h);
[I, ~] = find(h<0.5, 1, 'first');
fwhm = w(I)/pi;
fwhm = fwhm*param.fmax/(2*pi);
