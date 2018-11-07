%% Test Channel and Link classes
clc;
clear;
close all;

%% Default properties
channelArray = Channel('centerFrequency', -50e9, 'firFactor', 0.9, ...
    'powerdBm', -3);
% linkArray = Link('S', 0.06);
% linkArray = Link('gamma', 0);
linkArray = Link();

%% Other property values
channelArray(2) = Channel('modulation', '16QAM', 'symbolRate', 32e9, ...
    'centerFrequency', 0e9, 'firFactor', 0.2, ...
    'minSamplePerSymbol', 16, 'minNumberSymbol', 2^14, 'powerdBm', -3);

% linkArray(2) = Link('spanLength', 80e3, 'D', 1.6e-5, 'S', 0, ...
%     'alphadB', 0.25e-3, 'gamma', 1.5e-3, 'NFdB', 6, 'DCFLength', 80e3, ...
%     'numberSteps', 200);
% linkArray(2) = Link('S', 0.06);
% linkArray(2) = Link('gamma', 0);
linkArray(2) = Link();

%% Partial default values
channelArray(3) = Channel('centerFrequency', 50e9, ...
    'minSamplePerSymbol', 4, 'powerdB', -3);

% linkArray(3) = Link('DCFLength', 90e3);
% linkArray(3) = Link('S', 0.06);
% linkArray(3) = Link('gamma', 0);
linkArray(3) = Link();

%% Wrong values throws an error
% channelArray(4) = Channel('modulation', '32QAM');

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

% plot(t, abs(sp.txSignalTime))

%% Plot spectrum domain signal
% figure;
% plot(omega/2/pi, 10*log10(abs(sp.txSignalSpectrum).^2*1e12))

%% Simulate
sp.simulate();

%% Received spectrum
% figure;
% plot(omega/2/pi, 10*log10(abs(sp.currentSignalSpectrum).^2*1e12))
%% Plot spectrum domain signal after transmission
% figure;
% plot(omega/2/pi, abs(sp.channelArray(3).rxTime).^2)

%% Plot constellation diagram of received signals
figure; hold on;
plot(sp.channelArray(1).rxSymbolMatched, '.')
plot(sp.channelArray(1).rxCloudCenter, 'x')

figure; hold on;
plot(sp.channelArray(2).rxSymbolMatched, '.')
plot(sp.channelArray(2).rxCloudCenter, 'x')

figure; hold on;
plot(sp.channelArray(3).rxSymbolMatched, '.')
plot(sp.channelArray(3).rxCloudCenter, 'x')

%% Plot rotated constellation diagrams
figure; hold on;
plot(sp.channelArray(2).rxSymbolRotated, '.')
plot(sp.channelArray(2).txSymbolMatched, 'x')