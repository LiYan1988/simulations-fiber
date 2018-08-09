function xt = ift(xf, df)
% Inverse Fourier transform
xt = fftshift(fft(fftshift(xf)))*sqrt(df/(2*pi));