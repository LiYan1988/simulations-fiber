%% Prepare simulateSingleQAM on Glenn
clc;
close all;
clear;

%%
folderName = fullfile(pwd, 'cluster');
if ~exist(folderName, 'dir')
    mkdir(folderName)
end

copyfile('../../Channel.m', folderName)
copyfile('../../Link.m', folderName)
copyfile('../../SinglePolarization.m', folderName)
copyfile('simulateSingleQAM.m', folderName)
copyfile('run_sbatch.py', folderName)

%% Variables
% Split the variables into groups of 16, so that the resources in Glenn is
% used optimally
cpuPerNode = 16;
powerQAM = -20:1:10;
symbolRate = [32e9, 64e9];
powerArray = {};
symbolRateArray = {};
k = 1;
l = 1;
for n = 1:length(powerQAM)
    for m = 1:length(symbolRate)
        powerArray{k, l} = num2str(powerQAM(n));
        symbolRateArray{k, l} = num2str(symbolRate(m), '%d');
        l = l+1;
        if l>cpuPerNode
            l = 1;
            k = k+1;
        end
    end
end
tmp = strcat(' "', powerArray, ';', symbolRateArray, '" ');
variableArray = {};
for n = 1:size(tmp, 1)
    variableArray{n} = strcat(tmp{n, :});
end

%% Prepare bash files
for n = 1:size(powerArray, 1)
%     copyfile(fullfile(pwd, 'simulateSingleQAM.sh'), ...
%         fullfile(folderName, sprintf('simulateSingleQAM%d.sh', n)));
modifyBash(...
    fullfile(pwd, 'simulateSingleQAM.sh'), ...
    fullfile(folderName, sprintf('simulateSingleQAM%d.sh', n)), ...
    sprintf('simulateSingleQAM%d', n), ...
    [0, 2, 0, 0], ...
    variableArray{n});
end

%% Modify bash files
function modifyBash(fileNameSrc, fileNameDst, jobName, wallTime, variableArray)

fid = fopen(fileNameSrc, 'r');
tline = fgetl(fid);
i = 1;
A{i} = tline;
while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    A{i} = tline;
end
fclose(fid);

A{4} = sprintf('#SBATCH -J %s', jobName);
A{6} = sprintf('#SBATCH -t %1d-%02d:%02d:%02d', ...
    wallTime(1), wallTime(2), wallTime(3), wallTime(4));
A{7} = sprintf('#SBATCH -o %s.stdout', jobName);
A{8} = sprintf('#SBATCH -e %s.stderr', jobName);
A{16} = sprintf('array=( %s )', variableArray);
A{27} = sprintf('mkdir $SLURM_SUBMIT_DIR/%s', jobName);
A{28} = sprintf('cp -rf $TMPDIR/* $SLURM_SUBMIT_DIR/%s', jobName);

fid = fopen(fileNameDst, 'w');
for i = 1:numel(A)
    fprintf(fid, '%s\n', A{i});
    if A{i+1} == -1
        break
    end
end
fclose(fid);
end
