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

[pathNames, dirNames, fileNames] = dirwalk('resultsLevel1_100GHz');

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

%%
QAMSNRdB = mean(resultLevel1.SNRdB(:, [2, 4]), 2);
QAMSymbolRate = resultLevel1.symbolRate(:, 2);
QAMPowerdBm = resultLevel1.powerdBm(:, 2);
OOKSNRdB = min(resultLevel1.SNRdB(:, [1, 3, 5]), [], 2);

results = [QAMSymbolRate, QAMPowerdBm, QAMSNRdB, OOKSNRdB];
results = sortrows(results, [1, 2]);

meshYpower = reshape(results(:, 2), 31, []);
meshXsym = reshape(results(:, 1), 31, [])*1e-9;
meshQAM = reshape(results(:, 3), 31, []);
meshOOK = reshape(results(:, 4), 31, []);

%%
figure;
[M, c] = contour(meshXsym, meshYpower, meshQAM, [12, 11.5, 11, 10.5, 10]);
c.LineWidth = 2;
c.ShowText  = 'on';
c.LabelSpacing = 180;
xlabel('Symbol rate (GBaud)');
ylabel('Power (dBm)');
title('QAM SNR')
% savefig('QAMSNR.png')

%% 
figure;
[M, c] = contour(meshXsym, meshYpower, meshOOK, [-6, -2, 2, 6, 10, 14]);
c.LineWidth = 2;
c.ShowText  = 'on';
c.LabelSpacing = 180;
xlabel('Symbol rate (GBaud)');
ylabel('Power (dBm)');
title('OOK SNR')
% savefig('OOKSNR.png')
%%
figure;
meshTemp = meshQAM;
meshTemp(meshOOK<14) = 0;
[M, c] = contour(meshXsym, meshYpower, meshTemp, [12, 10]);
c.LineWidth = 2;
c.ShowText  = 'on';
c.LabelSpacing = 180;
xlabel('Symbol rate (GBaud)');
ylabel('Power (dBm)');
title('QAM SNR (OOK SNR > 14)')
% savefig('QAMSNR_OOKSNR>14.png')