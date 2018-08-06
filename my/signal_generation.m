clc;
clear;
close all;

% --- Frequency parameters
% 400 GHz frequency span, -fmax:fmax-df
fmax = 4e11; % [Hz]
% number of spectrum points
fn = 2^20;
% spectrum resolution
df = 2*fmax/fn; % [Hz]
% spectrum axis
f = df*[(0:fn/2-1) (-fn/2:-1)];
% ---

% --- Time parameters
% Time domain
% time window, -tmax:tmax-dt
tmax = pi/df; % [s]
% time resolution
dt = 2*tmax/fn; % [s]
% time axis
t = dt*(-fn/2:fn/2-1);
% ---

% --- Fiber parameters
% D = 17 [ps/nm/km]
% wavelength = 1.55 [um]
% speed_light = 2.99792458e5 [km/s]
% beta2 = -D*wavelength^2/(speed_light*2*pi) [ns^2/km]
distance = 100; % [km], span length
beta2 = -2.1683e-05; % [ns^2/km], GVD
gamma = 1.27e-3; % [(W*m)^-1], nonlinear coefficient of SMF 
alpha = 0.2; % [dB/km], attenuation in dB
alpha = log(10)*alpha/10; % [1/km] in linear
% ---

% --- Spatial parameters
Ld = dt^2/abs(beta2)*1e18; % [km], dispersion distance
zn = round(20*distance/Ld*gamma^2*1e6); % Number of z steps
dz = distance/zn; % [km], step size in z
% ---

% --- Generate a pulse
uu = sech(t*1e7).*exp(-0.5*1i*1*t.^2);
% ---

% plot(t, abs(uu).^2)

% --- Plot input pulse shape and spectrum
temp = fftshift(ifft(uu)).*(fn*dt)/sqrt(2*pi); % spectrum
figure; subplot(2,1,1);
plot(t, abs(uu).^2, '--k'); hold on;
axis([-tmax tmax 0 max(abs(uu).^2)*2]);
xlabel('Normalized Time');
ylabel('Normalized Power');
title('Input and Output Pulse Shape and Spectrum');
subplot(2,1,2);
plot(fftshift(f)/(2*pi), abs(temp).^2, '--k'); hold on;
axis([-fmax/1e5 fmax/1e5 0 max(abs(temp).^2)*2]);
xlabel('Normalized Frequency');
ylabel('Spectral Power');
% ---

% --- Store dispersive phase shifts to speedup code
dispersion = exp(0.5*1i*beta2*f.^2*dz); % phase factor
hhz = 1i*gamma^2*dz; % nonlinear phase factor
