%% Test Channel and Link classes
clc;
clear;
close all;

%% Default properties
channelArray = Channel('centerFrequency', -50e9);
linkArray = Link();

%% Other property values
channelArray(2) = Channel('modulation', '16QAM', 'symbolRate', 32e9, ...
    'centerFrequency', 0e9, 'firFactor', 0.2, ...
    'minSamplePerSymbol', 16, 'minNumberSymbol', 2^14);

linkArray(2) = Link('spanLength', 80e3, 'D', 1.6e-5, 'S', 0, ...
    'alphadB', 0.25e-3, 'gamma', 1.5e-3, 'NFdB', 6, 'DCFLength', 80e3, ...
    'numberSteps', 200);

%% Partial default values
channelArray(3) = Channel('centerFrequency', 50e9);

linkArray(3) = Link('DCFLength', 90e3);

%% Wrong values throws an error
% link(4) = Link('modulation', '32QAM');

%% SinglePolarization
sp = SinglePolarization('linkArray', linkArray, ...
    'channelArray', channelArray, 'randomSeed', 2938);
fclose(sp.logFid);
clearvars linkArray channelArray

%% Time and frequency axis
dt = sp.dt;
t = sp.t;
omega = sp.omega;
domega = sp.domega;
N = sp.N;

%% Plot dataTime 
% plot(t, abs(sp.channelArray(2).dataTime))
% figure;
% plot(t, abs(sp.channelArray(1).dataTime))
% figure;
% plot(t, abs(sp.channelArray(3).dataTime))
plot(t, abs(sp.txSignalTime))

%% Plot Spectrum domain signal
plot(omega/2/pi, 10*log10(abs(sp.txSignalSpectrum).^2*1e12))