function h = subfigure(M, N, k, fig)
% Resize a figure window to a certain fraction of the screen.
% Argument usage is similar to subplot.
%
% INPUTS:
%    M   : Vertical number of figures
%    N   : Horizontal number of figures
%    k   : Position of current window
%    fig : Figure number (optional)
% OUTPUTS:
%    h   : Figure handle
%
% Pontus Johannisson, 2009-04-28
% This software is distributed under the terms of the GNU General
% Public License version 2

if nargin < 3
    error('Usage: subfigure(M, N, k, fig)');
end

if nargin < 4
    h = figure;
else
    h = figure(fig);
end

if k > M*N
    error('Usage: k must be smaller than M*N');
end

row = rem(k - 1, M) + 1;
col = ceil(k/M);

x1 = (col - 1)/N;
x2 = (col - 0)/N;
y1 = 1 - (row - 1)/M;
y2 = 1 - (row - 0)/M;

units = get(h, 'units');
set(h, 'units', 'normalized', 'outerposition', [x1 y2 (x2 - x1) (y1 - y2)]);
axes; pause(0.01); % This is a workaround for a possible bug
                   % (sometimes the axes are not scaled properly
                   % otherwise.)
set(h, 'units', units);
