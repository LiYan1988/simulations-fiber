clc;
clear;
close all;

S = 1e6; % nubmer of symbols
M = 16;
symbols = qammod(randi([0, M-1], S, 1), M);

% Two ways of creating RRC
sps = 32;
span = 20;
roll_off = 0.1;

% use rcosdesign
rrc1 = rcosdesign(roll_off, span, sps);
txSig1 = upfirdn(symbols, rrc1, sps, 1);

% This method can remove overheads of filter and upsampling, good!
rrc2 = comm.RaisedCosineTransmitFilter('Shape','Square root', ...
    'RolloffFactor',roll_off, ...
    'FilterSpanInSymbols',span, ...
    'OutputSamplesPerSymbol',sps, 'Gain', sps);
txSig = rrc2(symbols);

%% Constellation diagram
% scatterplot(txSig(sps*(span/2)+1:end), sps)

%% Simple shift

%% 
rrc3 = comm.RaisedCosineReceiveFilter('InputSamplesPerSymbol',sps, ...
    'DecimationFactor',1, ...
    'Shape', 'Square root', ...
    'RolloffFactor', roll_off,...
    'FilterSpanInSymbols', span, 'Gain', 1/sps, ...
    'DecimationOffset', 0);
rxSig = rrc3(txSig);

a = reshape(rxSig, sps, []);
plot(a)

% delay = rrc3.FilterSpanInSymbols;

% Show symbols and rxSig are aligned
% plot(abs(rxSig(delay+1:end)))
% hold on
% plot(abs(symbols(1:end-delay)))

%% 
% figure;
% hold on;
% scatter(real(rxSig(delay+1:end)), imag(rxSig(delay+1:end)))
% scatter(real(symbols(1:end-delay)), imag(symbols(1:end-delay)))