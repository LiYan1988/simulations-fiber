%% Information retrieval from simulation results
clc;
close all;
clear;

%% Level 1, uncomment if want to copy all files into one folder
% 
% folderName = fullfile(pwd, 'postResultLevel1');
% if ~exist(folderName, 'dir')
%     mkdir(folderName)
% end

[pathNames, dirNames, fileNames] = dirwalk('resultsLevel1');

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
writetable(resultLevel1, 'resultsLevel1.csv')