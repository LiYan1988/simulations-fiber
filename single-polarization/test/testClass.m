%% Test Channel and Link classes
clc;
clear;
close all;

%% Default properties
channelArray = Channel('centerFrequency', -50e9, 'firFactor', 0.9, ...
    'powerdBm', -3);
linkArray = Link();

%% Other property values
channelArray(2) = Channel('modulation', '16QAM', 'symbolRate', 32e9, ...
    'centerFrequency', 0e9, 'firFactor', 0.2, ...
    'minSamplePerSymbol', 16, 'minNumberSymbol', 2^14, 'powerdBm', -3);
linkArray(2) = Link();

%% Partial default values
channelArray(3) = Channel('centerFrequency', 50e9, ...
    'minSamplePerSymbol', 4, 'powerdB', -3);
linkArray(3) = Link();

%% SinglePolarization
sp = SinglePolarization('linkArray', linkArray, ...
    'channelArray', channelArray, 'randomSeed', 2938);
clearvars linkArray channelArray

%% Simulate
sp.simulate();