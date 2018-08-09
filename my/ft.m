function xf = ft(xt, fn, df)
% Fourier transform

xf = fftshift(fft(fftshift(xt)))/fn/sqrt(df/(2*pi));
