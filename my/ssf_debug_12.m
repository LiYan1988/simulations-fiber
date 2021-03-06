clc;
clear;
close all;

% find a bug in OOk generation

%% Fiber Parameters
% -------------- Primary parameters
param.fmax = 2*pi*320*1e9; % [Hz], should be common multiples of all channels' bandwidths
param.fn = 2^16; % number of spectrum points

param.span_length = 82; % [km], span length
param.beta2 = -2.1683e-23; % [s^2/km], GVD, D=17 [ps/ns/km]
S = 0.06*1e6; % [s/(m^2*km)], third order dispersion
param.wavelength = 1550*1e-9; % [m], reference wavelength
param.light_speed = 2*1e8; % [m/s], speed of light in fiber
param.beta3 = (S-4*pi*param.light_speed/param.wavelength^3*param.beta2)*...
    (param.wavelength^2/(2*pi*param.light_speed))^2; % [s^3/km], third order dispersion
% d3=1i*beta3/6*(2*pi*FF).^3;
param.gamma = 1.27; % [(W*km)^-1], nonlinear coefficient of SMF
param.alpha = log(10)*0.2/10; % [1/km] in linear, 0.2 dB/km, positive number
param.dz = 0.2; % [km]
param.ase_exist = true;
param.nsp = 1.8; % [1] spontaneous emission factor, NF=5.5
param.h = 6.626*1e-34; % [J*s], [W*Hz^-2] Plank's constant
param.nu = param.light_speed/param.wavelength*1.5; % [Hz], light speed is in fiber, so 1.5 can bring it back to normal light speed

param.random_seed = 54790;

%% Test
power_dbm = 0;
param_mp = cell(length(power_dbm)); % [dBm], power of each channel
time_elapsed = zeros(size(power_dbm));

power_dbm_ook = -6;
power_dbm_qam = 0;
gauss_factor = 0.7;
bw_ghz_ook = 10*1e9;
bw_ghz_qam = 32*1e9;
grid_ghz = 200*1e9;

for k=1:length(power_dbm) % power of 16QAM
    fprintf('Iteration %d of %d started.\n', k, length(power_dbm))
    t = tic;
    for m=1:length(power_dbm) % power of OOK
        % Change channel uniform power
        param_temp = configure_channels_default_5(param, power_dbm_qam,...
            power_dbm_ook, bw_ghz_qam, bw_ghz_ook, grid_ghz);
        
        % Generate Signal
        param_temp = generate_signals(param_temp);
        
        % spectrum
        figure; hold on; grid on; box on;
        plot(param_temp.f_plot, abs(param_temp.data_mod_f_current).^2,...
            'linewidth', 2)
        set(gca, 'YScale', 'log')
        xlim([-250, 250])
        ylim([1e-20, 1e-11])
        xlabel('Frequency (GHz)')
        ylabel('PSD (W/Hz)')
        pbaspect([7 4 1])
        
%         figure; hold on; grid on; box on;
%         plot(param_temp.t_plot, param_temp.data_mod_t_channel{1})
        
        % eye diagram
        figure; hold on; grid on; box on;
        plot(mod(1:param_temp.fn, 4*param_temp.sample_per_symbol(1))*param_temp.dt*1e9, ...
            circshift(param_temp.data_mod_t_channel{1}, 32)/max(param_temp.data_mod_t_channel{1}(:)), '.')
        pbaspect([7, 4, 1])
        xlabel('Time (ns)')
        ylabel('Amplitude')
        title('OOK channel eye diagram')
        
        % Propagation through a link
        param_temp = simulate_link1(param_temp);
        
        param_mp{k, m} = param_temp;
    end
    time_elapsed(k) = toc(t);
    fprintf('Iteration %d of %d finished.\n', k, length(power_dbm))
    fprintf('This iteration takes %.2f minutes, total running time %.2f minutes.\n', ...
        time_elapsed(k)/60, sum(time_elapsed)/60)
    
    fprintf('------------------------------------------\n')
    fprintf('\n')
end

%% Save results
save('ssf_debug_12_variable_parameters.mat','-v7.3')

%% Plot results
scatterplot(param_temp.signal_received_constellation_derotate{1})

% figure; hold on; grid on; box on;
% plot(param_temp.f_plot, 10*log10(abs(param_temp.data_mod_f_in).^2))