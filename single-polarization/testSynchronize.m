%% Transmit a single QAM in 5 spans of fibers

clc;
close all;
clear;

%% Sweep power
spArray = SinglePolarization();
% QAMSnrdB = zeros(size(QAMPowerdBm));
for n=1:1
    linkArray = [Link('alphadB', 0.3e-3, 'gamma', 2e-3)];
    channelArray = [Channel(...
        'powerdBm', -14, ...
        'modulation', '16QAM', ...
        'firFactor', 0.2, ...
        'minNumberSymbol', 4096, ...
        'symbolRate', 29.99999e9); ...
%         Channel(...
%         'powerdBm', qamPowerdBm(n), ...
%         'modulation', '16QAM', ...
%         'firFactor', 0.2, ...
%         'minNumberSymbol', 4096, ...
%         'symbolRate', 30e9, ...
%         'centerFrequency', 50e9)...
        ];
    
    sp = SinglePolarization(...
        'simulationName', 'testSingleQAM', ...
        'simulationId', 1, ...
        'randomSeed', 0, ...
        'linkArray', linkArray, ...
        'channelArray', channelArray, ...
        'useParallel', false);
    sp.simulate();
    spArray(n) = sp;
    fprintf('SNR = %2.2f dB\n', sp.channelArray.SNRdB);
end

%% 
obj = spArray(1);
channelIdx = 1;
channel = obj.channelArray(1);
signal = channel.rxTime;
signal = [real(signal), imag(signal)];

q2 = synchronizeObj(obj, channel, signal, 0);

%%
function synchronize2(obj, channelIdx)
% Find the optimal place to sample the signal

channel = obj.channelArray(channelIdx);
% Copy received signal
rxSignal = channel.rxTime;
% Real and imaginary parts of the signal
signal = zeros(size(rxSignal, 1), 2);
signal(:, 1) = real(rxSignal);
signal(:, 2) = imag(rxSignal);

objfcn = @(x) -synchronizeObj(obj, channel, signal, x);

[xOpt,fval,exitflag,output,trials] = surrogateopt(___) 

channel.rxOptimalOffset

% Down sample received signal at the optimal offset
signal = downsample(signal, channel.actualSamplePerSymbol, ...
    channel.rxOptimalOffset-1);
channel.rxSymbol = signal(:, 1)+1i*signal(:, 2);
end

function q2 = synchronizeObj(obj, channel, signal, x)

    % 0<x<=channel.actualSamplePerSymbol
    if (x<=0) || (x>channel.actualSamplePerSymbol)
        q2 = -1;
        return
    end

    tmpSignal = downsample(signal, channel.actualSamplePerSymbol, ceil(x)-1);
    % Instruct kmeans to use parallel threads
    opts = statset('UseParallel', obj.useParallel);
    [~, tmpCenters, tmpDistances] = ...
        kmeans(tmpSignal, channel.constellationSize, ...
        'Display', 'off', 'maxiter', 1000, ...
        'Replicates', 4, 'Options', opts);
    
    tmpSignalPower = tmpCenters(:, 1)+1i*tmpCenters(:, 2);
    tmpSignalPower = abs(tmpSignalPower-tmpSignalPower.'); 
    tmpq = tmpSignalPower./(tmpDistances+tmpDistances.');
    q2 = mean(tmpq(:));
end