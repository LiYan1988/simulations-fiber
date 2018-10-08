function [U, f] = fourier_transform(u, t)
% Plot the amplitude and phase for the Fourier transform of the input
% signal using FFT and the proper scaling of the calculation.
%
% INPUTS:
%    u : The input signal [sqrt(W)]
%    t : Time vector [s]
% OUTPUTS:
%    none
%
% We assume that the time vector is monotonous and that the length
% is 2^N, where N is integer.
%
% Pontus Johannisson, 2009-05-08
% This software is distributed under the terms of the GNU General
% Public License version 2

% Calculate the period in time to get scaling of frequency vector
if nargin < 2
    dt = 1;
else
    dt = t(2) - t(1);
end
Nt = length(u);
T  = Nt*dt;
f  = fftshift(1/T*[0:Nt/2 - 1, -Nt/2:-1]);

% The Fourier transform amplitude should be scaled by dt. We also use
% fftshift to get the phase correct. (The time vector is assumed to be
% monotonous, not set up like frequency vector.)
U = fftshift(dt*fft(fftshift(u)));
