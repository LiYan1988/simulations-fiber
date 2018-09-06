function param = simulate_link1(param)

%% Propagation through a link
% 4 spans of 82 km SSMF followed by DCF with 80 km SSMF-equivalent
% compensation, plus one span of 82 km SSFM followed by DCF with 90 km
% SSMF-equivalent compensation

% FBG parameters
fbg_beta2 = -param.beta2;
fbg_beta3 = -param.beta3;
fbg_length_1 = 80; % 80 % [km]
fbg_length_2 = 90; % 90 % [km]

for n=1:4
    % Split Step Fourier in SSMF
    param = split_step_single_polarization(param);
    
    % FBG
    [signal_t_out, signal_f_out] = fbg_propagation(...
        param.data_mod_t_current, fbg_beta2, fbg_beta3, fbg_length_1, param);
    
    param.data_mod_t_current = signal_t_out;
    param.data_mod_f_current = signal_f_out;
end

% one span of 82 km SSFM, followed by DCF with 90 km SSMF-equivalent
% compression
param = split_step_single_polarization(param);
[signal_t_out, signal_f_out] = fbg_propagation(...
    param.data_mod_t_current, fbg_beta2, fbg_beta3, fbg_length_2, param);
param.data_mod_t_current = signal_t_out;
param.data_mod_f_current = signal_f_out;

%% Receiver 
% Recive signal and find centers of point clouds in constellations
dispersion_residual = param.span_length*5-fbg_length_1*4-fbg_length_2;
param = center_constellation(param, dispersion_residual);

% Rotate constellations back for visualization
param = derotate_constellation(param);

% Calculate SNR for each channel
param = calculate_snr(param);

% Calculate SER and EVM
param = calculate_ser_evm(param);
