function plot_iq(p, u, plot_opts)
% Plot the I and Q components of the input signal.
%
% INPUTS:
%    p         : The parameter struct
%    u         : The signal struct
%    plot_opts : Options to plot function
% OUTPUTS:
%   none
%
% Pontus Johannisson, 2010-01-22
% This software is distributed under the terms of the GNU General
% Public License version 2

if nargin < 3
    plot_opts = '';
end

ax(1) = subplot(2, 1, 1); plot(p.t_samp/1e-12, real(u), plot_opts);
ylabel('I channel [sqrt(W)]');
grid on;
box on;

ax(2) = subplot(2, 1, 2); plot(p.t_samp/1e-12, imag(u), plot_opts);
ylabel('Q channel [sqrt(W)]');
grid on;
box on;

xlabel('time [ps]');
linkaxes(ax, 'x');
zoom on;