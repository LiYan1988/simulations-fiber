clc;
clear;
close all;

% Write functions to implement:
% find centers of point clouds in the constellation diagram

%% Load data
load ssf_signal_polarization_5.mat

%% Test function
param = center_constellation(param);

%% Plot results
close all
for cidx= 1:param.channel_number%(param.channel_number+1)/2%1:param.channel_number
    figure;
    hold on;
    
    u = param.signal_received_constellation{cidx};
    plot(u(:, 1), u(:, 2), '.')
    
    v = param.cloud_centers{cidx};
    plot(v(:, 1), v(:, 2), 'x')
end

%% Save 
save ssf_signal_polarization_7.mat