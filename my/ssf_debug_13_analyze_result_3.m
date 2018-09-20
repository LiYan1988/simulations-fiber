clc;
clear;
close all;

% Analyze results from ssf_debug_13.m

%%
load ssf_debug_13_current_result.mat
% columns: ook power, qam power, grid_ghz, left ook snr, right ook snr, 
% qam snr, ook snr max, ook snr min, ook snr mean
results(:, 9) = (results(:, 7)+results(:, 8))/2;

%% 
% Remove some columns
% New columns: power ook, power qam, spacing GHz, snr qam, snr ook mean
results(:, [4, 5, 7, 8]) = [];
results(:, 3) = results(:, 3)/1e9;
results = array2table(results, 'VariableNames', {'power_OOK_dBm', ...
    'power_QAM_dBm', 'spacing_GHz', 'snr_QAM', 'snr_OOK_mean'});

writetable(results, 'snr_vs_power_spacing_results.csv')