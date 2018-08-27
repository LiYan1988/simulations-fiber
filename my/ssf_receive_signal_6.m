clc;
close all;
clear;

% Plot received signal after matching filter

%% Load data
load ssf_receive_signal_5.mat

%% Plot signals
cidx = 6;
x = param.data_mod_t_current;

% Downconverting
xt1 = x.*exp(1i*2*pi*param.center_frequency_channel(cidx).*param.t);
% Pass through matching filter
xt2 = upfirdn(xt1, param.filter_tx_channel{cidx}, 1, 1);
% Remove head and tail to make the constellation diagram clean
xt2 = xt2(param.filter_delay(cidx)+1:end-param.filter_delay(cidx));
% plot(param.f_plot, abs(ft(xt2, param.df)).^2)

scatterplot(...
    xt2(param.delay_filter_channel(cidx)+1:...
    end-param.delay_filter_channel(cidx)), ...
    param.sample_per_symbol(cidx), param.shift_channel_time(cidx))
