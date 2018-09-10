clc;
clear;
close all;

% Analyze results from ssf_debug_13.m

%%
load ssf_debug_13_current_result.mat
% columns: ook power, qam power, grid_ghz, ook1 snr, ook2 snr, qam snr,
% ook snr max, ook snr min
results(:, 9) = (results(:, 7)+results(:, 8))/2;


%%
figure; hold on; box on; grid on;
% plot(results(:, 7))
% plot(results(:, 8))
% plot(results(:, 6))

power_dbm_ook = -4;
power_dbm_qam = -4;
idx = (results(:, 1)==power_dbm_ook) & (results(:, 2)==power_dbm_qam);
a = results(idx, :);

plot(a(:, 3)/1e9, a(:, 6), 'displayname', '16 QAM')
% plot(a(:, 3), a(:, 7))
plot(a(:, 3)/1e9, a(:, 8), 'displayname', 'OOK')
legend('Location', 'Northwest')
title(sprintf('QAM power=%d dBm, OOK power=%d dBm', ...
    power_dbm_qam, power_dbm_ook))
xlabel('Channel spacing (GHz)')
ylabel('Linear SNR')