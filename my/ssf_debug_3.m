clc;
clear;
close all;

load matlab.mat

cidx=1;
xt = param.data_mod_t_current.*exp(1i*2*pi*param.center_frequency_channel(cidx).*param.t);
xf = ft(xt, param.df);
figure;
hold on;
plot(param.f_plot, 10*log10(abs(xf).^2))
xf = xf.*param.f_mask;
plot(param.f_plot, 10*log10(abs(xf).^2))
xt = ift(xf, param.df);

plot(param.f_plot, 10*log10(abs(xf).^2))