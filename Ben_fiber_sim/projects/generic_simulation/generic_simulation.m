function p = generic_simulation(p)
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
p = get_parameters(p);

% Generate and modulate light
u0 = light_signal(p);
measure_power(u0, 'laser output');

% Set up the electric driving signals
p.electric = p.electric1; e1 = electric_signal(p);
p.electric = p.electric2; e2 = electric_signal(p);

u0 = iq_modulator(p, u0, [e1; e2]);
measure_power(u0, 'after mod.  ');

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

% Propagate through a fiber
p.smf.delta_z = 1/max(max(abs(u0).^2))/p.smf.gamma/p.steps_per_L_NL;
p.fiber = p.smf;
u = nlse_spm(p, u0);
measure_power(u, 'after prop. ');

if isfield(p.fig, 'spectrum2')
    figure(p.fig.spectrum2);
    plot_spectrum(p, u);
    title('Spectrum after propagation');
end

if isfield(p.fig, 'eye2')
    figure(p.fig.eye2);
    plot_eye(p, u);
    title('Eye diagram after propagation');
end

if isfield(p.fig, 'const2')
    figure(p.fig.const2);
    plot_constellation(p, u);
    title('Constellation after propagation');
end
