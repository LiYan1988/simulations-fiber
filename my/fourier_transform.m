clc;
clear;
close all;

df = 2*pi*0.01; % Hz
fmax = pi*100; % Hz
f = -fmax:df:fmax-df;
fn = length(f);
dt = pi/fmax; % s
tmax = pi/df; % s
t = -tmax:dt:tmax-dt;
xt = cos(2*pi*10*t);
% Fourier transform
xf = fftshift(ifft(fftshift(xt)))/sqrt(df/(2*pi));
% Inverse Fourier transform
xt2 = fftshift(fft(fftshift(xf)))*sqrt(df/(2*pi));
xf2 = fftshift(ifft(fftshift(xt2)))/sqrt(df/(2*pi));
figure;
subplot(2, 1, 1)
plot(t, xt)
hold on;
plot(t, xt2, '--')
subplot(2, 1, 2)
plot(f/(2*pi), abs(xf))
hold on;
plot(f/(2*pi), abs(xf2))

Pt = norm(xt)^2/fn;
Pf = norm(xf)^2*df/(2*pi);
Pt2 = norm(xt2)^2/fn;
Pf2 = norm(xf2)^2*df/(2*pi);
disp(Pf)
disp(Pt)
disp(Pf2)
disp(Pt2)
%%
% M = 25;
% px = zeros(1, M);
% pxi = zeros(1, M);
% pxf = zeros(1, M);
% for n=1:M
%     N = 2^n;
%     x = randn(N, 1);
%     xf = fft(x)/sqrt(N);
%     xi = ifft(x)*sqrt(N);
%     
%     px(n) = norm(x);
%     pxf(n) = norm(xf);
%     pxi(n) = norm(xi);
% end
% 
% %%
% figure;
% hold on;
% plot(1:M, log2(pxf./px), 'o')
% plot(1:M, log2(px./pxi))
% % plot(1:M, 2.^(1:M), 'o')
