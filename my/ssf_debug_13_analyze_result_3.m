clc;
clear;
close all;

% Analyze results from ssf_debug_13.m

%%
load ssf_debug_13_current_result.mat
% columns: ook power, qam power, grid_ghz, left ook snr, right ook snr,
% qam snr, ook snr max, ook snr min, ook snr mean
results(:, 9) = (results(:, 7)+results(:, 8))/2;
results_ser(:, 4) = 0.5*(results_ser(:, 1)+results_ser(:, 3));
ser = array2table(results_ser(:, [2, 4]), 'VariableNames', ...
    {'ser_QAM', 'ser_OOK_mean'});
%%
% Remove some columns
% New columns: power ook, power qam, spacing GHz, snr qam, snr ook mean
results(:, [4, 5, 7, 8]) = [];
results(:, 3) = results(:, 3)/1e9;
results = array2table(results, 'VariableNames', {'power_OOK_dBm', ...
    'power_QAM_dBm', 'spacing_GHz', 'snr_QAM', 'snr_OOK_mean'});

writetable(results, 'snr_vs_power_spacing_results.csv')

%% What OOK power is suitable? -2 dBm
% Plot OOK SNR vs OOK power
snr_part = results(results{:, 'spacing_GHz'}==50, :);
ser_part = ser(results{:, 'spacing_GHz'}==50, :);

figure;
hold on;
grid on;
box on;
for power_qam = -4:1:3
    idx = snr_part{:, 'power_QAM_dBm'}==power_qam;
    plot(snr_part{idx, 'power_OOK_dBm'}, ...
        snr_part{idx, 'snr_OOK_mean'}, ...
        'displayname', sprintf('QAM power= %d dBm', power_qam))
end
legend()
plot(-10:1:3, 20*ones(14, 1), ...
    'linewidth', 2, 'linestyle', '--', 'color', 'r', ...
    'displayname', 'SNR threshold')
xlabel('OOK power (dBm)')
ylabel('OOK SNR')
% conclusion: OOK power = -2 dBm is optimal

%% What QAM power is suitable?
figure;
hold on;
grid on;
box on;
for power_ook = -6:1:-2
    idx = snr_part{:, 'power_OOK_dBm'}==power_ook;
    plot(snr_part{idx, 'power_QAM_dBm'}, ...
        snr_part{idx, 'snr_QAM'}, ...
        'displayname', sprintf('OOK power= %d dBm', power_ook))
end
plot(-10:1:3, 32.6*ones(14, 1), ...
    'linewidth', 2, 'linestyle', '--', 'color', 'r', ...
    'displayname', 'SNR threshold')
legend()
xlabel('QAM power (dBm)')
ylabel('QAM SNR')

% OOK should be less than -3 dBm to ensure QAM works well

%% Find feasible power range for both formats
idx = (snr_part{:, 'snr_OOK_mean'}>=25) & (snr_part{:, 'snr_QAM'}>=40);
snr_good = snr_part(idx, :);

%% Test the function
spacing_GHz = 40;
snr_QAM_th = 45;
snr_OOK_th = 25;
snr_good = find_good_snr(results, spacing_GHz, snr_QAM_th, snr_OOK_th);

%% Conclusion
% OOK power -7~-5, QAM power -2~0 @50GHz