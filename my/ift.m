function xt = ift(xf, fn, df)
% Inverse Fourier transform
xt = fftshift(ifft(fftshift(xf)))*fn*sqrt(df/(2*pi));