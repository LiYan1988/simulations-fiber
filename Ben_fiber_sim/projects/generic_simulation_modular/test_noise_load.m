function test_noise_load(snr)
% Test file
% Creates a signal, adds ASE noise and laser phase noise
% Saves the output file
% Optional: Process the noise-loaded output
%% Setting up paths
addpath('../../filter');
addpath('../../modulation');
addpath('../../nlse');
addpath('../../signal_plot');
addpath('../../transmission');
addpath('../../utils');

if nargin < 1
    snr = [0:20];
end

% Get parameters
p.timestamp = now();
p = get_parameters(p);

% Generate signal
[p, u0] = generate_wdm_signal(p);

% Save input
save(sprintf('%s%s_Parameters.mat', p.file.path, datestr(p.timestamp, 'yyyy-mm-dd_HH.MM.SS')), 'p')
save(sprintf('%s%s_Noise_load_input_waveform.mat', p.file.path, datestr(p.timestamp, 'yyyy-mm-dd_HH.MM.SS')), 'u0')

% Add TX laser phase noise
% phase_var = 2*pi*linewidth*p.dt_symb;
% u0 = add_phase_drift(u0, 2*pi*p.tx.linewidth*p.dt_symb);
% Add frequency offset
% u0 = u0.*exp(1j*p.t*0.2e9);

for i = 1:numel(snr)
    % Add ASE
    u = awgn(u0, snr(i), 'measured');

    % Add RX laser phase noise
%     u = add_phase_drift(u, 2*pi*p.tx.linewidth*p.dt_symb);

    % Save output
    save(sprintf('%s%s_Noise_load_output_waveform_%.2fdb.mat', p.file.path, datestr(p.timestamp, 'yyyy-mm-dd_HH.MM.SS'), snr(i)), 'u')

    % Coherent receiver
end

%% Cleaning up paths
rmpath('../../filter');
rmpath('../../modulation');
rmpath('../../nlse');
rmpath('../../signal_plot');
rmpath('../../transmission');
rmpath('../../utils');

end