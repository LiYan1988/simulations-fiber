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

[pathNames, dirNames, fileNames] = dirwalk('stderr');

errMsg = {};
resultLevel1 = struct();
i = 1;
for n = 1:length(fileNames)
    for m = 1:length(fileNames{n})
        errMsg{n, m} = fileread(fullfile(pathNames{n}, fileNames{n}{m}));
        fprintf('%s: %s', fileNames{n}{m}, errMsg{n, m})
    end
end

