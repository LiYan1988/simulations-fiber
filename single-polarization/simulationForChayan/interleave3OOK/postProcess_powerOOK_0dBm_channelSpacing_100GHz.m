%% Information retrieval from simulation results
% In Glenn command line, use 
% zip -r resultsLevel1.zip . -i *_resultsLevel1.mat
% to collect results
clc;
close all;
clear;

%% Level 1, uncomment if want to copy all files into one folder
% 
% folderName = fullfile(pwd, 'postResultLevel1');
% if ~exist(folderName, 'dir')
%     mkdir(folderName)
% end

channelSpacing = 100;
[pathNames, dirNames, fileNames] = dirwalk(sprintf('resultsLevel1_%dGHz', channelSpacing));

fileNameList = {};
resultLevel1 = struct();
i = 1;
for n = 1:length(fileNames)
    for m = 1:length(fileNames{n})
        fileNameList{i} = fullfile(pathNames{n}, fileNames{n}{m});
        %         copyfile(fileNameList{i}, folderName)
        tmp = load(fileNameList{i});
        if i==1
            resultLevel1 = tmp.result;
        else
            resultLevel1(i) = tmp.result;
        end
        i = i+1;
    end
end

resultLevel1 = struct2table(resultLevel1);
writetable(resultLevel1, sprintf('resultsLevel1_%dGHz.csv', channelSpacing))

%%
QAMSNRdB = mean(resultLevel1.SNRdB(:, [2, 4]), 2);
QAMSymbolRate = resultLevel1.symbolRate(:, 2);
QAMPowerdBm = resultLevel1.powerdBm(:, 2);
OOKSNRdB = min(resultLevel1.SNRdB(:, [1, 3, 5]), [], 2);

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
[M, c] = contour(meshXsym, meshYpower, meshQAM, [15, 14, 13, 12]);
c.LineWidth = 2;
c.ShowText  = 'on';
c.LabelSpacing = 180;
xlabel('QAM symbol rate (GBaud)');
ylabel('QAM power (dBm)');
title('QAM SNR (OOK power=0dBm)')

colormap(parula)
c = colorbar;
c.Label.String = 'SNR (dB)';
caxis([11, 16]) % [min-(max-min)/3, max+(max-min)/3]

savefig(fullfile(figureFolder, sprintf('QAMSNR_%dGHz.fig', channelSpacing)))
print(fullfile(figureFolder, sprintf('QAMSNR_%dGHz', channelSpacing)), '-dpng', '-r600')
print(fullfile(figureFolder, sprintf('QAMSNR_%dGHz', channelSpacing)), '-dpdf', '-r600', '-bestfit')

%% 
figure;
[M, c] = contour(meshXsym, meshYpower, meshOOK, [14, 16, 18, 20]);
c.LineWidth = 2;
c.ShowText  = 'on';
c.LabelSpacing = 180;
xlabel('QAM symbol rate (GBaud)');
ylabel('QAM power (dBm)');
title('OOK SNR (OOK power=0dBm)')

colormap(parula)
c = colorbar;
c.Label.String = 'SNR (dB)';
caxis([8, 22]) % [min-(max-min)/3, max+(max-min)/3]

savefig(fullfile(figureFolder, sprintf('OOKSNR_%dGHz.fig', channelSpacing)))
print(fullfile(figureFolder, sprintf('OOKSNR_%dGHz', channelSpacing)), '-dpng', '-r600')
print(fullfile(figureFolder, sprintf('OOKSNR_%dGHz', channelSpacing)), '-dpdf', '-r600', '-bestfit')

%%
OOKSNRth = 14;
figure;
meshTemp = meshQAM;
meshTemp(meshOOK<OOKSNRth) = 0;
[M, c] = contour(meshXsym, meshYpower, meshTemp, [15, 14, 13, 12]);
c.LineWidth = 2;
c.ShowText  = 'on';
c.LabelSpacing = 180;
xlabel('QAM symbol rate (GBaud)');
ylabel('QAM power (dBm)');
title(sprintf('QAM SNR (OOK power=0dBm, OOK SNR>%ddB)', OOKSNRth))

colormap(parula)
c = colorbar;
c.Label.String = 'SNR (dB)';
caxis([11, 16]) % [min-(max-min)/3, max+(max-min)/3]

savefig(fullfile(figureFolder, sprintf('QAMSNR_OOKSNR_greater_than_%ddB_%dGHz.fig', OOKSNRth, channelSpacing)))
print(fullfile(figureFolder, sprintf('QAMSNR_OOKSNR_greater_than_%ddB_%dGHz', OOKSNRth, channelSpacing) ), '-dpng', '-r600')
print(fullfile(figureFolder, sprintf('QAMSNR_OOKSNR_greater_than_%ddB_%dGHz', OOKSNRth, channelSpacing) ), '-dpdf', '-r600', '-bestfit')

