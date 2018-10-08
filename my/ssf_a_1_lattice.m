clc;
clear;
close all;

% OOK-QAM-OOK, power of OOK = -5 dBm

%% Fiber Parameters
% -------------- Primary parameters
param.fmax = 2*pi*1e9*(320); % [Hz], should be common multiples of all channels' bandwidths
param.fn = 2^17; % number of spectrum points

param.span_length = 82; % [km], span length
param.beta2 = -2.1683e-23; % [s^2/km], GVD, D=17 [ps/ns/km]
% S = 0.06*1e6; % [s/(m^2*km)], third order dispersion
param.wavelength = 1550*1e-9; % [m], reference wavelength
param.light_speed = 2.99792458*1e8; % [m/s], speed of light in vaccum
S = 4*pi*param.light_speed/param.wavelength^3*param.beta2;
param.beta3 = (S-4*pi*param.light_speed/param.wavelength^3*param.beta2)*...
    (param.wavelength^2/(2*pi*param.light_speed))^2; % [s^3/km], third order dispersion
% d3=1i*beta3/6*(2*pi*FF).^3;
param.gamma = 1.4; % [(W*km)^-1], nonlinear coefficient of SMF
param.alpha = log(10)*0.2/10; % [1/km] in linear, 0.2 dB/km, positive number
param.dz = 0.82; % [km]
param.ase_exist = true;
param.nsp = 1.2559; % [1] spontaneous emission factor, NF=5.563, nsp = 10^(NF/10)/2
param.h = 6.626*1e-34; % [J*s], [W*Hz^-2] Plank's constant
param.nu = param.light_speed/param.wavelength; % [Hz], light speed is in fiber, so 1.5 can bring it back to normal light speed

param.fbg_length_1 = 82; % [km], length of the first type DCF
param.fbg_length_2 = 82; % [km], length of the second type DCF

param.random_seed = 1; %54790;

%% Test
power_dbm_ook = -2;
power_dbm_qam = -10:1:6;
bw_hz_ook = 10*1e9; %(5:5:30)*1e9;
bw_hz_qam = (1:30)*1e9; %(5:5:50)*1e9;
grid_hz = 25*1e9;

% elapsed time
time_elapsed = zeros(length(power_dbm_qam), length(bw_hz_qam));

% cell to store results
param_mp = cell(length(power_dbm_qam), length(bw_hz_qam));

n = 1;
for k=1:length(power_dbm_qam)
    power_dbm_qam_temp = power_dbm_qam(k);
    for m=1:length(bw_hz_qam)
        fprintf('Iteration %d of %d started.\n', n, length(power_dbm_qam)*length(bw_hz_qam))
        fprintf('16QAM power %.1f dBm, baud rate %.1f G.\n', power_dbm_qam(k), bw_hz_qam(m)/1e9)
        t = tic;
        bw_hz_qam_temp = bw_hz_qam(m);
        param_temp = configure_channels_default_5(param, ...
            power_dbm_qam_temp, ...
            power_dbm_ook, ...
            bw_hz_qam_temp, bw_hz_ook, ...
            grid_hz);
        
        % Generate Signal
        param_temp = generate_signals(param_temp);
                plot(param_temp.f_plot, 10*log10(abs(param_temp.data_mod_f_in).^2*1e12))
        
        % Propagation through a link
        param_temp = simulate_link2(param_temp, 4);
        
        param_mp{k, m} = param_temp;
        time_elapsed(k, m) = toc(t);
        fprintf('Iteration finished.\n')
        fprintf('This iteration takes %.2f minutes, total running time %.2f minutes.\n', ...
            time_elapsed(k,m)/60, sum(time_elapsed(:))/60)
        fprintf('------------------------------------------\n')
        fprintf('\n')
        n = n+1;
        
        plot(param_temp.f_plot, abs(param_temp.data_mod_f_current).^2)
        plot(param_temp.t_plot, abs(param_temp.data_mod_t_current).^2)
        plot(param_temp.t_plot, abs(param_temp.data_mod_t_channel{2}).^2)
    end
end

save(sprintf('ssf_a_1_lattice-1qam2ook.mat'), '-v7.3')

%% Process
snr_results = zeros(510, 4);
n = 1;
for k=1:length(power_dbm_qam)
    for m=1:length(bw_hz_qam)
        snr_results(n, 1) = power_dbm_qam(k);
        snr_results(n, 2) = bw_hz_qam(m)/1e9;
        snr_results(n, 3) = param_mp{k, m}.snr_channel(2);
        snr_results(n, 4) = (param_mp{k, m}.snr_channel(1)+param_mp{k, m}.snr_channel(3))/2;
        n = n+1;
    end
end

save('ssf_a_1_lattice-1qam2ook.mat', 'snr_results', '-v7.3')