function simulateScenario(powerQAM1, powerQAM2, symbolRate1, symbolRate2, channelSpacing)
%% Simulation of 16QAM1+16QAM2
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
    Channel('modulation', '16QAM', ...
    'centerFrequency', 0e9, ...
    'symbolRate', symbolRate1, ...
    'powerdBm', powerQAM1, ...
    'minNumberSymbol', 2^14); ...
    Channel('modulation', '16QAM', ...
    'symbolRate', symbolRate2, ...
    'centerFrequency', channelSpacing, ...
    'powerdBm', powerQAM2, ...
    'minNumberSymbol', 2^14)];

%%
simulationName = sprintf('PQAM1_%d_PQAM2_%d_symbolRate1_%d_symbolRate2_%d_channelSpacing_%d', ...
    powerQAM1, powerQAM2, symbolRate1/1e9, symbolRate2/1e9, channelSpacing/1e9);

sp = SinglePolarization(...
    'simulationName', simulationName, ...
    'simulationId', 1, ...
    'linkArray', linkArray, ...
    'channelArray', channelArray, ...
    'useParallel', false);
sp.simulate();
sp.saveSimulationResult(true, false, false);

end