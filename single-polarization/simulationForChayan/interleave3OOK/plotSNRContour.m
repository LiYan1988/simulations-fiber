function plotSNRContour(resultLevel1, channelSpacing)
% Plot SNR contours, save figures as .fig, .png, and .pdf in a separate
% folder
% 
% INPUTs:
%   resultLevel1 is a table containing all results
% 
% OUTPUT:
%   no output
% 

%% Preprocess the result table
% QAMSNRdB = mean(resultLevel1.SNRdB(:, [2, 4]), 2);
% QAMSymbolRate = resultLevel1.symbolRate(:, 2);
% QAMPowerdBm = resultLevel1.powerdBm(:, 2);
% OOKSNRdB = min(resultLevel1.SNRdB(:, [1, 3, 5]), [], 2);
% 
% results = [QAMSymbolRate, QAMPowerdBm, QAMSNRdB, OOKSNRdB];
% results = sortrows(results, [1, 2]);

% nRow = length(unique(QAMPowerdBm));
% 
% meshYpower = reshape(results(:, 2), nRow, []);
% meshXsym = reshape(results(:, 1), nRow, [])*1e-9;
% meshQAM = reshape(results(:, 3), nRow, []);
% meshOOK = reshape(results(:, 4), nRow, []);


%% Create folder for figures
figureFolder = sprintf('figures_%dGHz', channelSpacing);
if ~exist(figureFolder, 'dir')
    mkdir(figureFolder)
end

%% SNR contour of QAM
figure;
[M, c] = contour(meshXsym, meshYpower, meshQAM, [15, 14, 13, 12]);
c.LineWidth = 2;
c.ShowText  = 'on';
c.LabelSpacing = 180;
xlabel('Symbol rate (GBaud)');
ylabel('Power (dBm)');
title('QAM SNR')
savefig(fullfile(figureFolder, sprintf('QAMSNR_%dGHz.fig', channelSpacing)))
print(fullfile(figureFolder, sprintf('QAMSNR_%dGHz', channelSpacing)), '-dpng', '-r600')
print(fullfile(figureFolder, sprintf('QAMSNR_%dGHz', channelSpacing)), '-dpdf', '-r600', '-bestfit')

