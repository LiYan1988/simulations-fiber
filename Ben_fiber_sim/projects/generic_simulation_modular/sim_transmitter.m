function p = sim_transmitter(slurm_id, p)
% A generic transmitter module. Creaes a parameter struct and generates a
% WDM signal. Should be paired with a '_link' function for actual
% transmission.
%
% INPUTS:
%    p : The parameter struct from caller
%
% Benjamin Foo, 2018-02-19
%
% Based on code by:
% Pontus Johannisson, 2010-04-19
% This software is distributed under the terms of the GNU General
% Public License version 2

%% Setting up paths
addpath('../../filter');
addpath('../../modulation');
addpath('../../nlse');
addpath('../../signal_plot');
addpath('../../transmission');
addpath('../../utils');

if ~exist('slurm_id', 'var')
    error('Specify the SLURM Job ID.')
end

%% Setting up parameter structure
if ~exist('p', 'var') || ~isfield(p, 'timestamp')
    p.timestamp = now(); % Get current time
end

p = get_parameters(p);

fprintf('%s, TX started: %s %s\n', datestr(p.timestamp, 'yyyy-mm-dd HH:MM:SS'), datestr(p.timestamp, 'yyyy-mm-dd_HH.MM.SS'), mfilename);

write_log(p, p.timestamp, sprintf('%s Started', mfilename));


%% Transmitter
[p, u0] = generate_wdm_signal(p);
p.link.wdm_pow=measure_power(u0);
save(sprintf('%s%s_Modulator_output_optical_waveform.mat', p.file.path, datestr(p.timestamp, 'yyyy-mm-dd_HH.MM.SS')), 'u0')

if p.flag.osa % Monitoring transmitted spectrum
    fig_handle = figure(80); clf(fig_handle); plot_spectrum(p, u0); title(fig_handle.CurrentAxes, 'TX Optical Spectrum');
    fig_handle.Position = [2560 850 700 500];
end

p = rmfield(p, {'f', 't'}); % Removing space-intensive fields. Re-calculate these when loading the parameter struct from file
save(sprintf('%s%s_Parameters.mat', p.file.path, datestr(p.timestamp, 'yyyy-mm-dd_HH.MM.SS')), 'p') % Saving parameter file

% Informing SLURM of the simulation ID (so that it knows what parameter
% files to load)
fid = fopen(sprintf('Job_%d_files.txt', slurm_id), 'w');
fprintf(fid, '%s%s', p.file.path, datestr(p.timestamp, 'yyyy-mm-dd_HH.MM.SS'));
fclose(fid);

fprintf('%s TX finished: %s %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), datestr(p.timestamp, 'yyyy-mm-dd_HH.MM.SS'), mfilename);
write_log(p, now, sprintf('%s Complete\n', mfilename));
if p.log.fid ~= -1 % If log is open, close it
    fclose(p.log.fid);
end

%% Removing paths
rmpath('../../filter');
rmpath('../../modulation');
rmpath('../../nlse');
rmpath('../../signal_plot');
rmpath('../../transmission');
rmpath('../../utils');
