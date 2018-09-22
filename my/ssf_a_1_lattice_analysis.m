clc;
clear;
close all;

% Analyze results, OOK-QAM-OOK, power of OOK = -5 dBm

%% Load data
load 'ssf_a_1_lattice_snr.mat'

power_qam_dbm = zeros(30, 17);
bw_qam_ghz = zeros(30, 17);
snr_qam = zeros(30, 17);
snr_ook = zeros(30, 17);
n = 1;
for k=1:17 % column, power
    for m=1:30 % row, baud rate
        power_qam_dbm(m, k) = snr_results(n, 1);
        bw_qam_ghz(m, k) = snr_results(n, 2);
        snr_qam(m, k) = snr_results(n, 3);
        snr_ook(m, k) = snr_results(n, 4);
        n = n+1;
    end
end

%% Plot contour
snr_qam2 = snr_qam;
snr_qam2(snr_ook<10) = 0;

figure; hold on;
grid on;
box on;
contour(bw_qam_ghz, power_qam_dbm, snr_qam2)
% contour(bw_qam_ghz, power_qam_dbm, snr_ook)
xlabel('Baud rate (GHz)')
ylabel('Power (dBm)')
