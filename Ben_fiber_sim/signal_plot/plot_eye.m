function plot_eye(p, u, plot_opts)
% Plot the eye diagram corresponding to the input field.
%
% INPUTS:
%    p         : The parameter struct
%    u         : The signal struct
%    plot_opts : Options to plot function
% OUTPUTS:
%   none
%
% Pontus Johannisson, 2009-05-06
% This software is distributed under the terms of the GNU General
% Public License version 2

if nargin < 3
    plot_opts = '';
end

% Set up the relative index vector for the current symbol slot and the
% neighboring two half slots.
rel_idx = (1 - p.samp_per_symb/2):(p.samp_per_symb + p.samp_per_symb/2);

% Corresponding time vector in ps
t_plot = (rel_idx - 1)*p.dt_samp/1e-12;

% Do second subplot first to make the title string occur at the top
% of the figure.

% Q channel
subplot(2, 1, 2); hold on;
for k = 1:p.N_symb
    u_idx = p.samp_per_symb*(k - 1) + rel_idx;
    if (k == 1) || (k == p.N_symb)
        u_idx = mod(u_idx - 1, p.N_samp) + 1;
    end
    plot(t_plot, imag(u(:, u_idx)), plot_opts);
end
xlabel('time [ps]')
ylabel('Q channel')
grid on;
zoom on;
box on;

% I channel
subplot(2, 1, 1); hold on;
for k = 1:p.N_symb
    u_idx = p.samp_per_symb*(k - 1) + rel_idx;
    if (k == 1) || (k == p.N_symb)
        u_idx = mod(u_idx - 1, p.N_samp) + 1;
    end
    plot(t_plot, real(u(:, u_idx)), plot_opts);
end
ylabel('I channel')
grid on;
zoom on;
box on;
