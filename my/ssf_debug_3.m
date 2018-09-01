clc;
clear
close all;

load matlab.mat

% Recive signal and find centers of point clouds in constellations
param = center_constellation(param);

% Rotate constellations back for visualization
param = derotate_constellation(param);

%%
close all
for n=1:param.channel_number
    tmp = param.signal_rx_complex{n};
    tmpC = param.cloud_centers{n};
    figure;
    hold on;
    plot(real(tmp), imag(tmp), '.')
    plot(tmpC(:, 1), tmpC(:, 2), 'x', 'linewidth', 2)
end


for n=1:param.channel_number
    tmp = param.signal_received_constellation_derotate{n};
    tmpC = param.cloud_centers_derotation{n};
    figure;
    hold on;
    plot(tmp(:, 1), tmp(:, 2), '.')
    plot(tmpC(:, 1), tmpC(:, 2), 'x', 'linewidth', 2)
end