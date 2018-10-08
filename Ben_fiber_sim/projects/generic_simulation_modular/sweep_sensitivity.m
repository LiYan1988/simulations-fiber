function q = sweep_sensitivity(psa_flag, sig_pwr_rx_in, data_file_in, path_in)
% MATLAB script to run a sensitivity sweep.
% 
% Either generate a new noiseless signal, or load a noiseless signal from a
% prevoius simulation. Then, add noise and measure performance for both PIA
% and PSA cases.
%
% Written by Benjamin Foo, 2018-02-19
% Photonics Lab, Chalmers University of Technology

%% Setting up paths
addpath('../../filter');
addpath('../../modulation');
addpath('../../nlse');
addpath('../../signal_plot');
addpath('../../transmission');
addpath('../../utils');

if ~exist('psa_flag', 'var')
    error('Please provide the first argument. Specify a ''1'' for PSA result or ''0'' for PIA result')
end

if ~exist('path_in', 'var')
    path_in = '';
end

%% Generating noiseless signal
% job_id = 1;
% p = psa_single_span_transmitter(job_id);
% 
% for sweep_ind = 1:length(disp_array)
%     p = psa_single_span_link(job_id, 'spans', 1, 'power', power, 'convergence', 1000);
% end

%% Load a noiseless signal
try
%     load('../../results/2018-04-25_17.13.13_Fiber_optical_waveform_-10_dBm_1000_nlse_steps_span_1.mat', 'u', 'u_idler')    
%     load('../../results/2018-04-25_17.13.13_Parameters.mat', 'p')

%     load('../../results/single_span_psa_2018-04-30/2018-04-30_13.59.32_Fiber_optical_waveform_-5_dBm_1000_nlse_steps_10.0_km_precomp_span_1.mat', 'u', 'u_idler')    
%     load('../../results/single_span_psa_2018-04-30/2018-04-30_13.59.32_Parameters.mat', 'p')    
    load([path_in, data_file_in], 'u', 'u_idler')
    
    path_idx = strfind(data_file_in, '/');
    if isempty(path_idx)
        tok_idx = strfind(data_file_in, '_');
        sim_id = data_file_in(1:tok_idx(2));
    else
        tok_idx = strfind(data_file_in(path_idx(end):end), '_');
        sim_id = data_file_in(1:path_idx(end)+tok_idx(2)-1);
    end
    load([path_in, sim_id, 'Parameters.mat'], 'p')
catch
    error('Could not find the requested files. Specify another simulation input.')
end
p.t = 0:p.dt_samp:(p.N_samp - 1)*p.dt_samp;
p.f = 1 /p.T_period*[(0:p.N_samp/2 - 1) (-p.N_samp/2:-1)];


%% Noise loading
sig_pwr_in = measure_power(u); % Signal power in dBm
try
    sig_pwr_rx = sig_pwr_rx_in; % Target received signal power
catch
    sig_pwr_rx = -40:-4:-52; 
end
vac_noise_pwr_dbm = 10*log10(p.const.h*(p.const.nu)*p.f_samp/2*1e3); % Vacuum noise power [dBm]

atten = sig_pwr_rx - sig_pwr_in;

evm_db = NaN(1, numel(atten));
gmi = NaN(1, numel(atten));
ber = NaN(1, numel(atten));

demux = struct('delta_f', p.rx.filter_bandwidth_hz, 'f0', 0, 'f', p.f); % Demultiplexing filter

for ind = 1:length(atten)
    u_sig_n = edfa(u, atten(ind));
%     u_sig_n = add_noise(u_sig_n, vac_noise_pwr_dbm); % Noise loading
    
    if psa_flag
        u_idl_n = edfa(u_idler, atten(ind));
%         u_idl_n = add_noise(u_idl_n, vac_noise_pwr_dbm); % Noise loading
    end
    
%     figure; plot_spectrum(p, u_sig_n);  % Testing noise loading function. Noise PSD should change proportionally to total noise power
%     figure; plot_spectrum(p, u_idl_n); pause; % Testing noise loading function. Noise PSD should change proportionally to total noise power

    if psa_flag
        % PSA pre-amplifier
%         % With reduced NF (possibly counts same improvement twice)
%         u_sig_n = edfa(u_sig_n, p.smf.L*p.smf.att_db, p.link.edfa_nf_db-3, p.const.nu, p.dt_samp);
%         u_idl_n = edfa(u_idl_n, p.smf.L*p.smf.att_db, p.link.edfa_nf_db-3, p.const.nu, p.dt_samp);       
%         [u_sig_n, ~] = coherent_superposition(u_sig_n, u_idl_n, 360);

        % With increased gain (developed independently; may be completely wrong, particularly in long-haul transmission)
        [u_sig_n, ~] = coherent_superposition(u_sig_n, u_idl_n, 360);        
        u_sig_n = 2*u_sig_n; % 6-dB PSA gain from coherent superposition
        u_sig_n = edfa(u_sig_n, p.smf.L*p.smf.att_db-10*log10(4), p.link.edfa_nf_db, p.const.nu, p.dt_samp);
        
    else
        % PIA pre-amplifier
        u_sig_n = edfa(u_sig_n, p.smf.L*p.smf.att_db, p.link.edfa_nf_db, p.const.nu, p.dt_samp); 
    end
    
    %% Processing results
    p.rx.chan = round(p.link.N_chan/2);
    u_rx = filter_ideal_bp(demux, u_sig_n.*exp(-1j*p.w_index.sig(p.rx.chan)*p.t));
    p.rx.edc_L = 0; % Set this to 0 for a PSA case, dumbass!
    p_rx(ind) = coherent_receiver(p, u_rx);
    evm_db(ind) = p_rx(ind).rx.evm_db;
    gmi(ind) = p_rx(ind).rx.gmi;
    ber(ind) = p_rx(ind).rx.ber;
    
end

q.evm_db = evm_db;
q.gmi = gmi;
q.ber = ber;
q.pwr_rx = sig_pwr_rx;

figure(200); clf; fig_evm = axes('NextPlot', 'add');
figure(201); clf; fig_gmi = axes('NextPlot', 'add');
figure(202); clf; fig_ber = axes('YScale', 'log');
plot(fig_evm, q.pwr_rx, evm_db, 'b.'); ylabel(fig_evm,'EVM (dB)'); xlabel(fig_evm,'Received signal power (dBm)')
plot(fig_gmi, q.pwr_rx, gmi, 'rx'); ylabel(fig_gmi,'GMI'); xlabel(fig_gmi,'Received signal power (dBm)'); ylim([0, 2.2]) % Limit set for QPSK
semilogy(fig_ber, q.pwr_rx, ber, 'bx--'); ylabel(fig_ber,'BER'); xlabel(fig_ber,'Received signal power (dBm)'); grid(fig_ber, 'on'); ylim([5e-5, 6e-1]) % Limit set to useful region

%% Cleaning up paths
rmpath('../../filter');
rmpath('../../modulation');
rmpath('../../nlse');
rmpath('../../signal_plot');
rmpath('../../transmission');
rmpath('../../utils');

end