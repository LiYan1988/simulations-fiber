%% Transmit a single OOK in 5 spans of fibers

clc;
close all;
clear;

%% Configure link and channel
ookPowerdBm = [-20:2:10];

linkArray = [Link(); Link(); Link(); Link(); Link()];
% channelArray = [Channel()];

%% Initialize and run simulation
% sp = SinglePolarization(...
%     'simulationName', 'testSingleOOK', ...
%     'simulationId', 1, ...
%     'randomSeed', 0, ...
%     'linkArray', linkArray, ...
%     'channelArray', channelArray);

% sp.simulate();

%% Sweep power
ookSnrdB = zeros(size(ookPowerdBm));
for n=1:length(ookPowerdBm)
%     fprintf('OOK power = %d dBm\n', ookPowerdBm(n));
    sp = SinglePolarization(...
    'simulationName', 'testSingleOOK', ...
    'simulationId', 1, ...
    'randomSeed', 0, ...
    'linkArray', linkArray, ...
    'channelArray', Channel('powerdBm', ookPowerdBm(n)));
    sp.simulate();
    ookSnrdB(n) = sp.channelArray(1).SNRdB;
    fprintf('OOK power = %d dBm, SNR = %2.2d dB\n', ookPowerdBm(n), ookSnrdB(n));
end

%%
figure;
hold on;
grid on;
box on;
plot(ookPowerdBm, ookSnrdB, 'o', 'linewidth', 2, 'displayname', 'SNR')
xlabel('Power (dBm)')
ylabel('SNR (dB)')
legend()