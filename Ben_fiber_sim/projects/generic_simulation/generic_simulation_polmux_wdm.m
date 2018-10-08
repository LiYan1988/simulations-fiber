function p = generic_simulation_polmux_wdm(p, varargin)
% Run a generic simulation. Modulate with an I/Q-modulator with
% unequally spaced electric driving signals. Propagate through a fiber
% and plot the final result.
%
% INPUTS:
%    p : The parameter struct from caller
%
% Benjamin Foo, 2018-02-15
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

%% Setting up parameter structure
if ~exist('p', 'var') || ~isfield(p, 'timestamp')
    p.timestamp = now(); % Get current time
end

p = get_parameters_polmux_wdm(p);

% Dynamic parameter allocation using keyval pairs in the function call
% Particularly useful for parameter sweeps
if nargin > 1
    for sweep_ind = 1:numel(varargin)
        try key = char(varargin(sweep_ind));
            switch lower(key)
                case 'power'
                    % Change PER CHANNEL launch power [dBm]
                    value = cell2mat(varargin(sweep_ind+1));
                    if isnumeric(value)
                        p.link.pwr_per_chan = value; 
                    end                    
                case 'convergence'
                    % Change number of steps per nonlinear length
                    value = cell2mat(varargin(sweep_ind+1));
                    if isnumeric(value)
                        p.steps_per_L_NL = value; 
                    end
                case 'spans'
                    % Change number of steps per nonlinear length
                    value = cell2mat(varargin(sweep_ind+1));
                    if isnumeric(value)
                        p.link.N_span = value; 
                        p.rx.edc_L = p.smf.L*p.link.N_span; % Remember to update any dependencies!
                    end
                otherwise
                    warning('Unknown command ''%s''.', key);
            end
        catch
            % This is the value part. Do nothing.
        end
    end
end

fprintf('%s, Simulation started: %s %s\n', datestr(p.timestamp, 'yyyy-mm-dd HH:MM:SS'), datestr(p.timestamp, 'yyyy-mm-dd_HH.MM.SS'), mfilename);
write_log(p, p.timestamp, sprintf('%s Started', mfilename));
write_log(p, now, sprintf('Setup: Steps per nonlinear length %d ', p.steps_per_L_NL)); % For convergence testing


%% Transmitter
[p, u0] = generate_wdm_signal(p);
save(['../results/', datestr(p.timestamp, 'yyyy-mm-dd_HH.MM.SS'), '_Modulator_output_optical_waveform.mat'], 'u0');

if p.flag.osa % Monitoring transmitted spectrum
    fig_handle = figure(80); clf(fig_handle); plot_spectrum(p, u0); title(fig_handle.CurrentAxes, 'Input Optical Spectrum');
    fig_handle.Position = [2060 550 560 420];
end

if p.flag.b2b || (p.link.N_span == 0)
    %% Coherent receiver
    demux = struct('delta_f', p.rx.filter_bandwidth_hz, 'f0', 0, 'f', p.f);
    p.rx.chan = 2;
    if p.rx.chan > p.link.N_chan
        p.rx.chan = p.link.N_chan;
    end
    
    u_rx = filter_ideal_bp(demux, u0.*exp(-1j*p.w_index.sig(p.rx.chan)*p.t));
    p = coherent_receiver(p, u_rx);
    return
end

if isfield(p, 'fig') && isfield(p.fig, 'spectrum') && p.flag.osa
    figure(p.fig.spectrum);
    plot_spectrum(p, u0);
    title('Spectrum after modulation');
end

if isfield(p, 'fig') && isfield(p.fig, 'eye')
    figure(p.fig.eye);
    plot_eye(p, u0);
    title('Eye diagram after modulation');
end

if isfield(p, 'fig') && isfield(p.fig, 'const')
    figure(p.fig.const);
    plot_constellation(p, u0);
    title('Constellation after modulation');
end

%% Propagate through a fiber
% Set powers
u = set_power(u0, p.link.pwr_per_chan + 10*log10(p.link.N_chan));
measure_power(u, 'into link.  ')
write_log(p, now, sprintf('Target per-channel launch power is %d dBm', p.link.pwr_per_chan));

% Transmit through fiber spans
for span_ind = 1:p.link.N_span
    rng(mod(p.rng.seed + span_ind, 2^32-1)); % Re-setting RNG every span
    p.smf.delta_z = 1/max(max(abs(u).^2))/p.smf.gamma/p.steps_per_L_NL;
    p.fiber = p.smf;
    if isfield(p.flag, 'gpu') && p.flag.gpu
        u = nlse_spm_gpu(p, u);
    else
        u = nlse_spm(p, u);
    end
    measure_power(u, sprintf('after span %d', span_ind));

    u = edfa(u, p.smf.L*p.smf.att_db, p.link.edfa_nf_db, p.const.c/p.const.lambda, p.dt_samp);
    
    % Saving fiber output
    save(['../results/', datestr(p.timestamp, 'yyyy-mm-dd_HH.MM.SS'), '_Fiber_optical_waveform_', num2str(span_ind*p.smf.L/1e3), '_km.mat'], 'u')
    write_log(p, now, sprintf('SMF: Propagated %d km', span_ind*p.smf.L/1e3));    
end

write_log(p, now, sprintf('Transmission through %d km SMF completed', p.link.N_span*p.smf.L/1e3));

if isfield(p, 'fig') && isfield(p.fig, 'spectrum2') && p.flag.osa
    figure(p.fig.spectrum2);
    plot_spectrum(p, u);
    title('Spectrum after propagation');
end

if isfield(p, 'fig') && isfield(p.fig, 'eye2')
    figure(p.fig.eye2);
    plot_eye(p, u);
    title('Eye diagram after propagation');
end

if isfield(p, 'fig') && isfield(p.fig, 'const2')
    figure(p.fig.const2);
    plot_constellation(p, u);
    title('Constellation after propagation');
end

%% Coherent receiver
demux = struct('delta_f', p.rx.filter_bandwidth_hz, 'f0', 0, 'f', p.f);
p.rx.chan = 5;
if p.rx.chan > p.link.N_chan
    p.rx.chan = p.link.N_chan;
end

u_rx = filter_ideal_bp(demux, u.*exp(-1j*p.w_index.sig(p.rx.chan)*p.t));
p = coherent_receiver(p, u_rx);
p = rmfield(p, {'fiber', 'f', 't'}); % Removing useless fields
save(sprintf('..\\results\\%s_Parameters.mat', datestr(p.timestamp, 'yyyy-mm-dd_HH.MM.SS')), 'p');

fprintf('%s Simulation finished: %s %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), datestr(p.timestamp, 'yyyy-mm-dd_HH.MM.SS'), mfilename);
write_log(p, now, sprintf('%s Complete\n', mfilename));
if p.log.fid ~= -1 % If log is open, close it
    fclose(p.log.fid);
end

%% Removing paths
rmpath('../filter');
rmpath('../modulation');
rmpath('../nlse');
rmpath('../signal_plot');
rmpath('../transmission');
rmpath('../utils');
