clc;
clear;
close all;

% Write functions to implement:
% de-rotate the constellation disgram back for visualization purpose

%% Load data
load ssf_signal_polarization_7.mat

%% Test function
param = derotate_constellation(param);