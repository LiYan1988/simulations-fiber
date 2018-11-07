%% Search for the optimal factor of Gauss FIR
% The optimal bt is 1.
clc;
clear;
close all;

%% Define links and channels
% S = 0.06 ps/(nm^2*km) = 60 s/m^3
linkArray = [...
    Link('spanLength', 82e3, 'DCFLength', 80e3, 'S', 60); ...
    Link('spanLength', 82e3, 'DCFLength', 80e3, 'S', 60); ...
    Link('spanLength', 82e3, 'DCFLength', 80e3, 'S', 60); ...
    Link('spanLength', 82e3, 'DCFLength', 80e3, 'S', 60); ...
    Link('spanLength', 82e3, 'DCFLength', 90e3, 'S', 60)];

% channelArray = [Channel('modulation', 'OOK')];

%% Simulate
firFactor = 0.5:0.05:1;
ookPower = -10:1:10;
spArray(length(firFactor), length(ookPower)) = SinglePolarization();
for n = 1:length(firFactor)
    parfor m = 1:length(ookPower)
        simulationId = (n-1)*length(firFactor)+m;
        ookChannel = Channel('firFactor', firFactor(n), ...
            'powerdBm', ookPower(m));
        spArray(n, m) = SinglePolarization('channelArray', ookChannel, ...
            'linkArray', linkArray, ...
            'simulationName', 'testOptimalGaussFIR', ...
            'simulationId', simulationId);
        spArray(n, m).simulate();
        fprintf('bt = %.2f, power = %d dBm, SNR = %2.2f dB\n', ...
            firFactor(n), ookPower(m), spArray(n, m).channelArray.SNRdB);
    end
end

%%
SNRdBm = zeros(size(spArray));
for n = 1:length(firFactor)
    for m = 1:length(ookPower)
        SNRdBm(n, m) = spArray(n, m).channelArray.SNRdB;
    end
end

SNRdBmMax = max(SNRdBm, [], 2);
plot(firFactor, SNRdBmMax)

%% Eye 
plotEye(spArray(5, 10), 1, 4, 'rx', 1)