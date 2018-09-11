clc;
clear;
close all;

% Analyze results from ssf_debug_13.m

%% Test
power_dbm_ook = -10:1:6;
power_dbm_qam = -10:1:6;
bw_ghz_ook = 10*1e9; %(5:5:30)*1e9;
bw_ghz_qam = 32*1e9; %(5:5:50)*1e9;
grid_ghz = (30:5:200)*1e9;

% columns: ook power, qam power, grid_ghz, ook1 snr, ook2 snr, qam snr
results = zeros(0, 6);
results_ser = zeros(0, 3);

%%
u = 0;
v = 0;
listing = dir;
for n=1:length(listing)
    if length(listing(n).name)>17 && strcmp(listing(n).name(1:17), 'ssf_debug_13_qam_')
        t = tic;
        fprintf('%s\n', listing(n).name)
        
        load(listing(n).name)
        for g=1:length(param_mp)
            u = u+1;
            results(u, 1) = power_dbm_ook_temp;
            results(u, 2) = power_dbm_qam_temp;
            results(u, 3) = param_mp{g}.spectrum_grid_size;
            results(u, 4) = param_mp{g}.snr_channel(1);
            results(u, 5) = param_mp{g}.snr_channel(3);
            results(u, 6) = param_mp{g}.snr_channel(2);
            results_ser(u, :) = param_mp{g}.ser_channel;
        end
    end
end

%%
results = sortrows(results, [1, 2, 3]);
results(:, 7) = max(results(:, 4), results(:, 5));
results(:, 8) = min(results(:, 4), results(:, 5));
figure; hold on; box on; grid on;
plot(results(:, 7))
plot(results(:, 8))
plot(results(:, 6))

save ssf_debug_13_current_result.mat results results_ser