function [x] = dispersive_propagation(x, omega_vec, DL)
% Apply chromatic dispersion to the input signal

% INPUTS:
%    x         : Input signal
%    omega_vec : Angular frequency vector
%    DL        : Dispersion parameter x Distance
% OUTPUTS:
%    x : Output signal
%
% The linear Schrödinger equation describes dispersive propagation.
% This can be described as an all-pass filter, as described, e.g.,
% by Agrawal's book.
%
% The DL input argument is the product of the dispersion parameter and
% the propagation distance. We keep all parameters in SI-unit without
% prefix. This means that using the dispersion of a typical SMF D = 16
% ps/(nm km) = 16e-6 s/m^2 and a propagation distance L = 1000 km =
% 1e6 m gives DL = 16e-6*1e6 = 16.
%
% Copyright 2011 Pontus Johannisson
% This software is distributed under the terms of the GNU General
% Public License version 2

c          = 2.99792458e8;     % Speed of light [m/s]
lambda     = 1.55e-6;          % Wavelength [m]
D_to_beta2 = -lambda^2/2/pi/c; % * D [s/m^2]
beta2_L    = D_to_beta2*DL;
for k = 1:size(x, 1);
    x(k, :) = ifft(fft(x(k, :)).*exp(1i*1/2*beta2_L*omega_vec.^2));
end;
