%% Transmit a single OOK in 5 spans of fibers

clc;
close all;
clear;

%% Configure link and channel
qamPowerdBm = [-10:2:10];

linkArray = [Link(); Link(); Link(); Link(); Link()];

%% Sweep power
spArray = SinglePolarization();
% ookSnrdB = zeros(size(ookPowerdBm));
for n=5:length(qamPowerdBm)
%     fprintf('OOK power = %d dBm\n', ookPowerdBm(n));
    sp = SinglePolarization(...
    'simulationName', 'testSingleOOK', ...
    'simulationId', 1, ...
    'randomSeed', 0, ...
    'linkArray', linkArray, ...
    'channelArray', Channel(...
        'powerdBm', qamPowerdBm(n), ...
        'modulation', '16QAM', ...
        'firFactor', 0.2, ...
        'minNumberSymbol', 4096, ...
        'symbolRate', 32e9));
    sp.simulate();
    spArray(n) = sp;
    fprintf('QAM power = %d dBm, SNR = %2.2f dB\n', qamPowerdBm(n), sp.channelArray.SNRdB);
end

%%
EVMArray = ones(1, 11);
SNRArray = ones(1, 11);
for n=1:11
    EVMArray(n) = -20*log10(spArray(n).channelArray.EVM);
    SNRArray(n) = spArray(n).channelArray.SNRdB;
end

figure;
hold on;
grid on;
box on;
plot(qamPowerdBm, EVMArray, 'o', 'linewidth', 2, 'displayname', 'EVM')
plot(qamPowerdBm, SNRArray, 'x', 'linewidth', 2, 'displayname', 'SNR')
xlabel('Power (dBm)')
ylabel('EVM (dB)')
legend()
