%% Test Channel and Link classes
clc;
clear;

%% Default properties
channelArray = Channel();
linkArray = Link();

%% Other property values
channelArray(2) = Channel('modulation', '16QAM', 'symbolRate', 32e9, ...
    'centerFrequency', 50e9, 'filterFactor', 0.2, ...
    'minSamplePerSymbol', 16);

linkArray(2) = Link('spanLength', 80e3, 'D', 1.6e-5, 'S', 0, ...
    'alphadB', 0.25e-3, 'gamma', 1.5e-3, 'NFdB', 6, 'DCFLength', 80e3, ...
    'numberSteps', 200);

%% Partial default values
channelArray(3) = Channel('centerFrequency', -50e9);

linkArray(3) = Link('DCFLength', 90e3);

%% Wrong values throws an error
% link(4) = Link('modulation', '32QAM');

%% SinglePolarization
sp = SinglePolarization('linkArray', linkArray, 'channelArray', channelArray);
fclose(sp.logFid);

%% Channels
disp(sp.fmax/1e9)