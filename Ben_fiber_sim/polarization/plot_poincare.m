function plot_poincare()
% Draw a Poincaré sphere, points on the sphere are plotted with
% plot_stokes.m
%
% INPUTS:
%    none
%
% Pontus Johannisson, 2011-05-12

N_theta  = 7;    % Number of theta values (latitudes)
N_phi    = 13;   % Number of phi values (longitudes)
N_points = 500;  % Number of points for the sphere circles
sc_col   = 0.4*[1 1 1]; % Sphere circle color
ax_l     = 1.3;  % Axis length
ax_w     = 1.4;  % Axis width
tp       = 1.15; % Axis marking position parameter
fs       = 14;   % Fontsize for axis markings

% Draw latitudes
theta = linspace(0,   pi, N_theta);
phi   = linspace(0, 2*pi, N_points);
for k = 1:length(theta);
    plot3(sin(theta(k)).*cos(phi), ...
          sin(theta(k)).*sin(phi), ...
          cos(theta(k)).*ones(size(phi)), ...
          'Color', sc_col);
    hold on;
end;

% Draw longitudes
theta = linspace(0,   pi, N_points);
phi   = linspace(0, 2*pi, N_phi);
for k = 1:length(phi);
    plot3(sin(theta).*cos(phi(k)), ...
          sin(theta).*sin(phi(k)), ...
          cos(theta), ...
          'Color', sc_col);
end;

% Draw the axes
plot3([-ax_l ax_l], [0 0], [0 0], 'k', 'LineWidth', ax_w);
plot3([0 0], [-ax_l ax_l], [0 0], 'k', 'LineWidth', ax_w);
plot3([0 0], [0 0], [-ax_l ax_l], 'k', 'LineWidth', ax_w);

% Add labels
text(tp*ax_l, 0, 0, 'S_1', 'FontSize', fs);
text(0, tp*ax_l, 0, 'S_2', 'FontSize', fs);
text(0, 0, tp*ax_l, 'S_3', 'FontSize', fs);

axis equal;
axis off;
view(135, 30); % Change to a better view angle
rotate3d on;