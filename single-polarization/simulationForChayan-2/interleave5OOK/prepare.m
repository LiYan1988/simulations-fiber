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
copyfile('simulateScenario.m', folderName)
copyfile('run_sbatch.py', folderName)

%% Variables
% Split the variables into groups of 16, so that the resources in Glenn is
% used optimally
cpuPerNode = 16;
powerQAM = -20:1:10;
powerOOK = [0];
% symbolRate = [32e9, 64e9];
% channelSpacing = [50e9, 100e9, 150e9, 200e9];
symbolRate = (1:1:40)*1e9;
channelSpacing = [50e9];
wallTime = [0, 10, 0, 0];

parameterArray = combvec(powerQAM, powerOOK, symbolRate, channelSpacing);

% array=( "one;1;a;r" "two;2;b;s" "three;3;c;t" )
% powerQAM, powerOOK, symbolRate, channelSpacing
vstr = cell(ceil(length(parameterArray)/cpuPerNode), cpuPerNode);
rowIdx = 1;
colIdx = 1;
for n = 1:length(parameterArray)
    vstr{rowIdx, colIdx} = sprintf(' "%d;%d;%d;%d" ', ...
        parameterArray(1, n), parameterArray(2, n), ...
        parameterArray(3, n), parameterArray(4, n));
    colIdx = colIdx+1;
    if colIdx>cpuPerNode
        rowIdx = rowIdx+1;
        colIdx = 1;
    end
end

variableArray = cell(size(vstr, 1), 1);
for n=1:length(variableArray)
    variableArray{n} = strcat(vstr{n, :});
end

%% Prepare bash files
for n = 1:length(variableArray)
    modifyBash(...
        fullfile(pwd, 'simulateScenario.sh'), ...
        fullfile(folderName, sprintf('simulateScenario%d.sh', n)), ...
        sprintf('simulateScenario%d', n), ...
        wallTime, ...
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
A{28} = sprintf('cp -rf $TMPDIR/results/* $SLURM_SUBMIT_DIR/%s', jobName);

fid = fopen(fileNameDst, 'w');
for i = 1:numel(A)
    fprintf(fid, '%s\n', A{i});
    if A{i+1} == -1
        break
    end
end
fclose(fid);
end
