%% Test simulateSingleQAM

clc;
close all;
clear;

%% 
folderName = 'simulateSingleQAM';
if ~exist(folderName, 'dir')
    mkdir(folderName)
end

copyfile('Channel.m', folderName)
copyfile('Link.m', folderName)
copyfile('SinglePolarization.m', folderName)
copyfile('simulateSingleQAM.m', folderName)
copyfile('simulateSingleQAM.sh', folderName)
