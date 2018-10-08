function plot_stokes(S, varargin)
% Plot a number of Stokes vectors on the Poincaré sphere, which is
% drawn with plot_poincare.m
%
% INPUTS:
%    S        : The Stokes vectors (real, 3xN)
%    varargin : Specification like to plot.m (optional)
% OUTPUT:
%    none
%
% Pontus Johannisson, 2009-10-14

% Plot the points
if nargin > 1;
    plot3(S(1, :), S(2, :), S(3, :), varargin{:});
else
    plot3(S(1, :), S(2, :), S(3, :), 'k.', 'MarkerSize', 20);
end;
axis equal;
axis off;
view(135, 30); % Change to a better view angle
rotate3d on;

return
% Test cases, using Damask, p. 35
close all;
plot_poincare();
plot_stokes(jones2stokes([1;  0]),         'b*', 'MarkerSize', 10);
plot_stokes(jones2stokes([0;  1]),         'bo', 'MarkerSize', 10);
plot_stokes(jones2stokes([1;  1]/sqrt(2)), 'g*', 'MarkerSize', 10);
plot_stokes(jones2stokes([1; -1]/sqrt(2)), 'go', 'MarkerSize', 10);
plot_stokes(jones2stokes([1;  i]/sqrt(2)), 'r*', 'MarkerSize', 10);
plot_stokes(jones2stokes([1; -i]/sqrt(2)), 'ro', 'MarkerSize', 10);
plot_stokes(jones2stokes([cos(pi/8); sin(pi/8)]));
