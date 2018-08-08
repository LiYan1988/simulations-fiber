clc;
clear;
close all;

N = 2^20;
x = randn(N, 1); % x is a signal
y = ifft(x); 
z = fft(x);
px = norm(x)^2/N; % power of x in time domain
py = norm(y)^2; % power of x in frequency domain
pz = norm(z)^2/N^2; % I don't know...

% Two things to note:
% 1. ifft(x) and fft(x) are not unitary
% 2. ifft(x) is more suitable to be used as a Fourier transform, because
% the results can be directly used to plot the spectral power density
% figure. If use fft(x), I should normalize with a factor of 1/N.
% 3. In split step Fourier, we use
% sum(abs(x).^2)/N==sum(abs(x_f).^2)*domega/(2*pi) to conserve power in
% both domains. 
% 
% Fourier transform
% u_f = fftshift(ifft(u_t))*sqrt(df/(2*pi)) = fftshift(ifft(u_t))/sqrt(fn*dt)