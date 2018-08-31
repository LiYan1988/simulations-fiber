clc;
close all;
clear;

%% Simulate OOK
% 1. Understand OOK transmitter and receiver
% 2. No optics

%% Signal Generation
% yout = upfirdn(xin,h,p,q) specifies the integer downsampling factor q, 
% where q has a default value of 1. The length of the output, yout, is 
% ceil(((length(xin)-1)*p+length(h))/q)
rng(24)
bits_length = 100;
% bits = randi([0, 1], bits_length, 1);
bits = [1, zeros(1, bits_length-2), 1];

samples_per_symbol = 100;
symbols_in_filter = 6;
h = gaussdesign(1, symbols_in_filter, samples_per_symbol);

x0 = upfirdn(bits, [1], samples_per_symbol);
x = upfirdn(bits, h, samples_per_symbol);
x = x/max(x);

% Filter length is N, N is an odd number.
% Samples per symbol is M, M is an even number.
% The group delay of the filter is (N-1)/2
% But we want to keep the first half symbol, so we remove 1:(N-1)/2-M/2
% And start the signal from (N-1)/2-M/2+1
% The signal length should be multiplied by M, so the rest are removed as 
% tail
s = ((symbols_in_filter/2-0.5)*samples_per_symbol+1):...
    ((symbols_in_filter/2-0.5)*samples_per_symbol+bits_length*samples_per_symbol);
x1 = x(s);

hold on;
plot(x0)
plot(x)
plot(x1)