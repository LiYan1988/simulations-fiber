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

channelSpacing = 50;
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

%% Remove duplicate rows
resultsNew = results(1, :);
for n=2:length(results)
    tmp = results(n, :);
    if (tmp(1)==resultsNew(end, 1)) && (tmp(2)==resultsNew(end, 2))
        continue;
    else
        resultsNew(end+1, :) = tmp;
    end
end

save('resultsNew_50GHz.mat', 'resultsNew')

%% Extract columns
meshYpower = reshape(resultsNew(:, 2), 31, []);
meshXsym = reshape(resultsNew(:, 1), 31, [])*1e-9;
meshQAM = reshape(resultsNew(:, 3), 31, []);
meshOOK = reshape(resultsNew(:, 4), 31, []);

%%
figureFolder = sprintf('figures_%dGHz', channelSpacing);
if ~exist(figureFolder, 'dir')
    mkdir(figureFolder)
end

%%
% fig = figure;
% [M, c] = contour(meshXsym, meshYpower, meshQAM, [13, 12, 10, 11]);
% c.LineWidth = 2;
% c.ShowText  = 'on';
% c.LabelSpacing = 180;
% clabel(M,c,'Interpreter','latex')
% xlabel('QAM symbol rate (GBaud)', 'Interpreter', 'latex', 'fontsize', 12);
% ylabel('QAM power (dBm)', 'Interpreter', 'latex', 'fontsize', 12);
% grid on;
% pbaspect([7 4 1])
% 
% fig.PaperPosition = [0 0 7 4];
% 
% savefig(fullfile(figureFolder, sprintf('QAMSNR_%dGHz.fig', channelSpacing)))
% print(fullfile(figureFolder, sprintf('QAMSNR_%dGHz', channelSpacing)), '-dpng')
% print(fullfile(figureFolder, sprintf('QAMSNR_%dGHz', channelSpacing)), '-dpdf', '-r600', '-fillpage')

%% 
% fig = figure;
% [M, c] = contour(meshXsym, meshYpower, meshOOK, [6, 8, 10, 12, 14]);
% c.LineWidth = 2;
% c.ShowText  = 'on';
% clabel(M,c,'Interpreter','latex')
% c.LabelSpacing = 180;
% xlabel('QAM symbol rate (GBaud)', 'Interpreter', 'latex', 'fontsize', 12);
% ylabel('QAM power (dBm)', 'Interpreter', 'latex', 'fontsize', 12);
% grid on;
% pbaspect([7 4 1])
% fig.PaperPosition = [0 0 7 4];
% 
% savefig(fullfile(figureFolder, sprintf('OOKSNR_%dGHz.fig', channelSpacing)))
% print(fullfile(figureFolder, sprintf('OOKSNR_%dGHz', channelSpacing)), '-dpng', '-r600')
% print(fullfile(figureFolder, sprintf('OOKSNR_%dGHz', channelSpacing)), '-dpdf', '-r600', '-fillpage')

%%
OOKSNRth = 14;
figure;
meshTemp = meshQAM;
meshTemp(meshOOK<OOKSNRth) = 0;
[M, c] = contour(meshXsym, meshYpower, meshTemp, [13, 12, 11, 10]);
c.LineWidth = 2;
c.ShowText  = 'on';
c.LabelSpacing = 330;
clabel(M,c,'Interpreter','latex');
xlabel('QAM symbol rate (GBaud)', 'Interpreter','latex', 'fontsize', 12);
ylabel('QAM power (dBm)', 'Interpreter','latex', 'fontsize', 12);
grid on;
pbaspect([7 4 1]);
fig.PaperPosition = [0 0 7 4];

savefig(fullfile(figureFolder, sprintf('QAMSNR_OOKSNR_greater_than_%ddB_%dGHz.fig', OOKSNRth, channelSpacing)))
print(fullfile(figureFolder, sprintf('QAMSNR_OOKSNR_greater_than_%ddB_%dGHz', OOKSNRth, channelSpacing) ), '-dpng', '-r600')
print(fullfile(figureFolder, sprintf('QAMSNR_OOKSNR_greater_than_%ddB_%dGHz', OOKSNRth, channelSpacing) ), '-dpdf', '-r600', '-fillpage')

