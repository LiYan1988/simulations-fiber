function p = plot_constellation(p, u, line_color, samp_color, ms)
% Plot the constellation diagram corresponding to the input fields.
%
% INPUTS:
%    p          : The parameter struct
%    u          : The signal (complex)
%    line_color : Line color (1 x 3, optional)
%    samp_color : Sample color (1 x 3, optional)
%    ms         : Marker size (integer, optional)
% OUTPUTS:
%   none
%
% Use the input arguments line_color and samp_color to get
% user-defined colors. They are both (1 x 3) matrices.
%
% Pontus Johannisson, 2009-04-28
% This software is distributed under the terms of the GNU General
% Public License version 2

% Adding possibility to mark the first sample (for spotting "edge effects")
p.mark_first_sample = 0;

if nargin < 3 % Default line color
    line_color = [0, 0, 1];
end

if nargin < 4 % Default center sample color
    samp_color = [1, 0, 0];
end

if nargin < 5 % Default marker size
    ms = 10;
end

for k = 1:size(u, 1)
    if p.mark_first_sample
        plot(real(u(k, 1)),   imag(u(k, 1)),   'gx');
        plot(real(u(k, end)), imag(u(k, end)), 'g+');
    end
    hold on;
    plot(real(u(k, :)), imag(u(k, :)), 'Color', line_color);
    % The sampling instant is in the middle of the two center samples
    plot(real(u(k, p.idx_symb)), imag(u(k, p.idx_symb)), '.', ...
         'Color', samp_color, 'MarkerSize', ms);
    axis equal;
    grid on;
    zoom on;
    box on;
end
