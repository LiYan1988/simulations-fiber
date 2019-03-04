%% This file:
% Data of 200 GHz simulation is lost, or maybe in another disk. But the
% final result csv is here. In this file, use the final csv file to plot figures.

%% Information retrieval from simulation results
% In Glenn command line, use 
% zip -r resultsLevel1.zip . -i *_resultsLevel1.mat
% to collect results
clc;
close all;
clear;
channelSpacing = 200;

%% Read results
resultLevel1 = readtable('resultsLevel1_200GHz.csv');

%%
QAMSNRdB = mean([resultLevel1.SNRdB_2, resultLevel1.SNRdB_4], 2);
QAMSymbolRate = resultLevel1.symbolRate_2;
QAMPowerdBm = resultLevel1.powerdBm_2;
OOKSNRdB = min([resultLevel1.SNRdB_1, resultLevel1.SNRdB_3, resultLevel1.SNRdB_5], [], 2);

results = [QAMSymbolRate, QAMPowerdBm, QAMSNRdB, OOKSNRdB];
results = sortrows(results, [1, 2]);

nRow = length(unique(QAMPowerdBm));

meshYpower = reshape(results(:, 2), nRow, []);
meshXsym = reshape(results(:, 1), nRow, [])*1e-9;
meshQAM = reshape(results(:, 3), nRow, []);
meshOOK = reshape(results(:, 4), nRow, []);

%%
figureFolder = sprintf('figures_%dGHz', channelSpacing);
if ~exist(figureFolder, 'dir')
    mkdir(figureFolder)
end

%%
figure;
[M, c] = contour(meshXsym, meshYpower, meshQAM, [17, 15, 9]);
c.LineWidth = 2;
c.ShowText  = 'on';
c.LabelSpacing = 180;
xlabel('QAM symbol rate (GBaud)');
ylabel('QAM power (dBm)');
title('QAM SNR (OOK power=0dBm)')

colormap(parula)
c = colorbar;
c.Label.String = 'SNR (dB)';
caxis([6, 20]) % [min-(max-min)/3, max+(max-min)/3]

savefig(fullfile(figureFolder, sprintf('QAMSNR_%dGHz.fig', channelSpacing)))
print(fullfile(figureFolder, sprintf('QAMSNR_%dGHz', channelSpacing)), '-dpng', '-r600')
print(fullfile(figureFolder, sprintf('QAMSNR_%dGHz', channelSpacing)), '-dpdf', '-r600', '-bestfit')

%% 
figure;
[M, c] = contour(meshXsym, meshYpower, meshOOK, [10, 14, 16, 18]);
c.LineWidth = 2;
c.ShowText  = 'on';
c.LabelSpacing = 180;
xlabel('QAM symbol rate (GBaud)');
ylabel('QAM power (dBm)');
title('OOK SNR (OOK power=0dBm)')

colormap(parula)
c = colorbar;
c.Label.String = 'SNR (dB)';
caxis([7, 23]) % [min-(max-min)/3, max+(max-min)/3]

savefig(fullfile(figureFolder, sprintf('OOKSNR_%dGHz.fig', channelSpacing)))
print(fullfile(figureFolder, sprintf('OOKSNR_%dGHz', channelSpacing)), '-dpng', '-r600')
print(fullfile(figureFolder, sprintf('OOKSNR_%dGHz', channelSpacing)), '-dpdf', '-r600', '-bestfit')

%%
OOKSNRth = 14;
figure;
meshTemp = meshQAM;
meshTemp(meshOOK<OOKSNRth) = 0;
[M, c] = contour(meshXsym, meshYpower, meshTemp, [17, 15, 9]);
c.LineWidth = 2;
c.ShowText  = 'on';
c.LabelSpacing = 180;
xlabel('QAM symbol rate (GBaud)');
ylabel('QAM power (dBm)');
title(sprintf('QAM SNR (OOK power=0dBm, OOK SNR>%ddB)', OOKSNRth))

colormap(parula)
c = colorbar;
c.Label.String = 'SNR (dB)';
caxis([6, 20]) % [min-(max-min)/3, max+(max-min)/3]

savefig(fullfile(figureFolder, sprintf('QAMSNR_OOKSNR_greater_than_%ddB_%dGHz.fig', OOKSNRth, channelSpacing)))
print(fullfile(figureFolder, sprintf('QAMSNR_OOKSNR_greater_than_%ddB_%dGHz', OOKSNRth, channelSpacing) ), '-dpng', '-r600')
print(fullfile(figureFolder, sprintf('QAMSNR_OOKSNR_greater_than_%ddB_%dGHz', OOKSNRth, channelSpacing) ), '-dpdf', '-r600', '-bestfit')

