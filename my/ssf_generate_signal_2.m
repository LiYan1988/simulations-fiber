clc;
clear;
close all;

%% Parameters

% Primary parameters
params.fmax = 2*pi*400*1e9; % [Hz]
params.fn = 2^22; % number of spectrum points

params.distance = 100; % [km], span length
params.beta2 = -2.1683e-05; % [ns^2/km], GVD
params.gamma = 1.27e-3; % [(W*m)^-1], nonlinear coefficient of SMF 
params.alpha = log(10)*0.2/10; % [1/km] in linear, 0.2 dB/km
params.zn = 10000; % number of steps per span

% Channel specific parameters, n channels should have n sets of parameters
params.bandwidth_channel = 32*1e9; % [Hz]

% Secondary parameters
params.df = 2*fmax/fn; % [Hz], frequency resolution
params.f = df*[(0:fn/2-1) (-fn/2:-1)]';

% Time domain parameters
params.tmax = pi/df; % [s]
params.dt = 2*tmax/fn; 
params.t = dt*(-fn/2:fn/2-1)';

params.dz = distance/zn; % [km], step size in z

params.samples_per_symbol = ceil(fmax/pi/params.bandwidth_channel);
params.symbol_number = floor(fn/params.samples_per_symbol);

