function simulateScenario(powerQAM, powerOOK, symbolRate, channelSpacing)
%% Simulation of OOK+16QAM+OOK
% Variables:
% QAM power: -20:1:10
% channel spacing: [50, 100, 150, 200] GHz
% QAM symbol rate: 1:1:channelSpacing+10 GHz
%
% Constant:
% OOK power: 0
% OOK symbol rate: 10 GHz

% Note: for OOK channel, its optical filter bandwidth is always 50 GHz

numberChannel = 9;

%% Define links and channels
% S = 0.06 ps/(nm^2*km) = 60 s/m^3
linkArray = [...
    Link('spanLength', 82e3, 'DCFLength', 80e3, 'S', 60); ...
    Link('spanLength', 82e3, 'DCFLength', 80e3, 'S', 60); ...
    Link('spanLength', 82e3, 'DCFLength', 80e3, 'S', 60); ...
    Link('spanLength', 82e3, 'DCFLength', 80e3, 'S', 60); ...
    Link('spanLength', 82e3, 'DCFLength', 90e3, 'S', 60)];

channelArray = repmat(Channel(), numberChannel, 1);
for n = 1:numberChannel
    if mod(n, 2)==1
        channelArray(n) = ...
            Channel('modulation', 'OOK', ...
            'centerFrequency', (n-(numberChannel+1)/2)*channelSpacing, ...
            'powerdBm', powerOOK, ...
            'minNumberSymbol', 2^12, ...
            'opticalFilterBandwidth', 50e9, ... % Note: for OOK channel, its optical filter bandwidth is always 50 GHz
            'firFactor', 1 ... % For OOK, firFactor=1 means a normal Gaussian pulse
            );
    else
        channelArray(n) = ...
            Channel('modulation', '16QAM', ...
            'centerFrequency', (n-(numberChannel+1)/2)*channelSpacing, ...
            'symbolRate', symbolRate, ...
            'powerdBm', powerQAM, ...
            'minNumberSymbol', 2^12, ...
            'opticalFilterBandwidth', channelSpacing, ... % for QAM, not necessary to specify optical filter bandwidth
            'firFactor', 0.2... % for QAM, firFactor is the roll-off factor which is 0.2
            );
    end
end

%%
simulationName = sprintf('PQAM_%d_POOK_%d_symbolRate_%d_channelSpacing_%d', ...
    powerQAM, powerOOK, symbolRate/1e9, channelSpacing/1e9);

simulationNumber = 4;
i = 1;
while(i<=simulationNumber)
    sp = SinglePolarization(...
        'simulationName', simulationName, ...
        'simulationId', i, ...
        'linkArray', linkArray, ...
        'channelArray', channelArray, ...
        'useParallel', false, ...
        'saveObject', false, ...
        'randomSeed', i);
    sp.simulate();
    [resultTemp, ~, ~] = sp.saveSimulationResult(false, false, false);
    if i==1
        resultArray = resultTemp;
    else
        resultArray(end+1) = resultTemp;
    end
    i = i+1;
end

resultAverage = resultArray(1);
resultAverage.SNR = mean(cat(3, resultArray(:).SNR), 3);
resultAverage.SNRdB = 10*log10(resultAverage.SNR);
resultAverage.EVM = mean(cat(3, resultArray(:).EVM), 3);
resultAverage.SER = mean(cat(3, resultArray(:).SER), 3);
resultAverage.BER = mean(cat(3, resultArray(:).BER), 3);
resultAverage.achievableDataRate = mean(cat(3, resultArray(:).achievableDataRate), 3);
resultAverage.runningTime = sum(cat(3, resultArray(:).runningTime), 3);

matFileLevel1 = sprintf('%s_%d_resultsLevel1.mat', sp.simulationName, 0);
matFileLevel1 = fullfile(sp.resultFolder, matFileLevel1);
save(matFileLevel1, 'resultAverage');

end