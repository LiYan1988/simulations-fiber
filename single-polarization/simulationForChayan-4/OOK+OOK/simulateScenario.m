function simulateScenario(powerOOK1, powerOOK2, symbolRate, channelSpacing)
%% Simulation of OOK1+OOK2
% Variables:
% QAM power: -20:1:10 
% OOK power: -20:1:10 
% channel spacing: [50, 100, 150, 200] GHz
% OOK symbol rate: 10 GHz

%% Define links and channels
% S = 0.06 ps/(nm^2*km) = 60 s/m^3
linkArray = [...
    Link('spanLength', 82e3, 'DCFLength', 80e3, 'S', 60); ...
    Link('spanLength', 82e3, 'DCFLength', 80e3, 'S', 60); ...
    Link('spanLength', 82e3, 'DCFLength', 80e3, 'S', 60); ...
    Link('spanLength', 82e3, 'DCFLength', 80e3, 'S', 60); ...
    Link('spanLength', 82e3, 'DCFLength', 90e3, 'S', 60)];

channelArray = [...
    Channel('modulation', 'OOK', ...
    'centerFrequency', 0e9, ...
    'symbolRate', symbolRate, ...
    'powerdBm', powerOOK1, ...
    'minNumberSymbol', 2^12, ...
    'opticalFilterBandwidth', 50e9, ...
    'firFactor', 1); ...
    Channel('modulation', 'OOK', ...
    'centerFrequency', channelSpacing, ...
    'symbolRate', symbolRate, ...
    'powerdBm', powerOOK2, ...
    'minNumberSymbol', 2^12, ...
    'opticalFilterBandwidth', 50e9, ...
    'firFactor', 1)];

%%
simulationName = sprintf('POOK1_%d_POOK2_%d_symbolRate_%d_channelSpacing_%d', ...
    powerOOK1, powerOOK2, symbolRate/1e9, channelSpacing/1e9);

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