function xf = ft(xt, df)
% Fourier transform

xf = fftshift(ifft(fftshift(xt)))/sqrt(df/(2*pi));
