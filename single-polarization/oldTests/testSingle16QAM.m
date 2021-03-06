%% Transmit a single QAM in 5 spans of fibers

clc;
close all;
clear;

%% Configure link and channel
qamPowerdBm = [-10:2:10];


%% Sweep power
spArray = SinglePolarization();
% QAMSnrdB = zeros(size(QAMPowerdBm));
for n=1:length(qamPowerdBm)
    linkArray = [Link(); Link(); Link(); Link(); Link()];
    channelArray = [Channel(...
        'powerdBm', qamPowerdBm(n), ...
        'modulation', '16QAM', ...
        'firFactor', 0.2, ...
        'minNumberSymbol', 4096, ...
        'symbolRate', 30e9); ...
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
        'useParallel', true);
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

%% Check Channel properties
% mc = ?Channel;
% propertyList = mc.PropertyList;
% for n=1:length(propertyList)
%     if strcmp(propertyList(n).Name, 'constellationSize')
%         disp(n)
%     end
% end