clc;
clear;
close all;

% Modify downconvert function

%% Load data
load ssf_receive_signal_5.mat

%% Plot signals
cidx = 6;
x = param.data_mod_t_current;

[xt, xp, padding_length] = downconvert(x, param, cidx);

%%
scatterplot(...
    xt(param.delay_filter_channel(cidx)+1:...
    end-param.delay_filter_channel(cidx)), ...
    param.sample_per_symbol(cidx), param.shift_channel_time(cidx))