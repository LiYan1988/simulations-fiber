function p = generic_simulation_polmux_full_out(p)
% Run a generic simulation. Modulate with an I/Q-modulator with
% unequally spaced electric driving signals. Propagate through a fiber
% and plot the final result.
%
% INPUTS:
%    p : The parameter struct from caller
%
% Pontus Johannisson, 2010-04-19
% This software is distributed under the terms of the GNU General
% Public License version 2

addpath('../../filter');
addpath('../../modulation');
addpath('../../nlse');
addpath('../../signal_plot');
addpath('../../transmission');
addpath('../../utils');

p.exists = []; % Make sure p exists

% Set up the parameter structure
p = get_parameters_polmux(p);

% Generate and modulate light
u0 = light_signal(p);
measure_power(u0, 'laser output');

% Set up the electric driving signals
p.electric = p.electric1; e1 = electric_signal(p);
p.electric = p.electric2; e2 = electric_signal(p);
p.electric = p.electric3; e3 = electric_signal(p);
p.electric = p.electric4; e4 = electric_signal(p);

u0_x = iq_modulator(p, u0/sqrt(2), [e1; e2]); % u0/sqrt(2) because light from the same laser is split into two paths for polmux operation (power is halved, field is divided by sqrt(2))
u0_y = iq_modulator(p, u0/sqrt(2), [e3; e4]);
u0=[u0_x; u0_y];
measure_power([u0; u0(1, :)+u0(2, :)], 'after mod.  ');

P_total = measure_power(u0(1, :)+u0(2, :));

if isfield(p.fig, 'spectrum')
    figure(p.fig.spectrum);
    plot_spectrum(p, u0);
    title('Spectrum after modulation');
end

if isfield(p.fig, 'eye')
    figure(p.fig.eye);
    plot_eye(p, u0);
    title('Eye diagram after modulation');
end

if isfield(p.fig, 'const')
    figure(p.fig.const);
    plot_constellation(p, u0);
    title('Constellation after modulation');
end

u0 = edfa(u0, p.signal.launch_power_dBm-P_total); % Noiseless EDFA to set launch power into link
measure_power([u0; u0(1, :)+u0(2, :)], 'into fiber. ');

% Propagate through a fiber
p.smf.delta_z = 1/max(max(abs(u0(1, :)).^2+abs(u0(2, :)).^2))/p.smf.gamma/p.steps_per_L_NL;
p.fiber = p.smf;
% u = nlse_spm(p, u0);
u = nlse_manakov_tod_all_steps(p, u0);
measure_power([u(end-2:end-1, :); u(end-2, :)+u(end-1, :)], 'after prop. ');

if isfield(p.fig, 'spectrum2')
    figure(p.fig.spectrum2);
    plot_spectrum(p, u(end-2:end-1, :));
    title('Spectrum after propagation');
end

if isfield(p.fig, 'eye2')
    figure(p.fig.eye2);
    plot_eye(p, u(end-2:end-1, :));
    title('Eye diagram after propagation');
end

if isfield(p.fig, 'const2')
    figure(p.fig.const2);
    plot_constellation(p, u(end-2:end-1, :));
    title('Constellation after propagation');
end

% Saving results for later use
current_dir = pwd; % Current directory
parameters = rmfield(p, 'fig'); % Removing figure handles from structure
parameters.signal.e_field_in = u0; % Adding the input electric field of the signal
parameters.signal.e_field_out = u; % Adding the output electric field of the signal

% Adding step size
if isinf(p.fiber.delta_z);
    parameters.nlse.steps = 1;
else
    parameters.nlse.steps = ceil(p.fiber.L/p.fiber.delta_z); % Number of iterations
end;

% Format for file with entire data structure: [Timestamp (year-month-day hour.minute.second)]_[Signal]_[Link]_[Power].mat
save_file_name = ['[', datestr(now, 'yyyy-mm-dd_HH.MM.SS'), ']_', num2str(p.N_symb), 'symb_', num2str(p.f_symb/1e9), 'Gbaud_', p.modulation.format,'_', num2str(p.smf.L/1e3), 'km_', num2str(p.signal.launch_power_dBm, '%.1f'), 'dBm.mat'];
% Check to make sure output directory exists
if ~exist('../../results', 'dir')
    mkdir('../../results')
end

save(['../../results/', save_file_name], 'parameters')
display(['Saved output to: ', save_file_name])