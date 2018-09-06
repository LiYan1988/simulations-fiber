clc;
clear;
% close all;

% Calculate EVM and BER

%%
load debug_9_2.mat

plot(param.f_plot, 10*log10(abs(param.data_mod_f_current).^2))