function [x, phi] = add_phase_drift(x, phase_var)
% Introduce a phase drift modeled as a Wiener process
%
% INPUTS:
%    x         : Input vector
%    phase_var : Variance for Wiener process
% OUTPUTS:
%    x   : Output vector
%    phi : The phase [rad]
%
% Copyright 2011 Pontus Johannisson
% This software is distributed under the terms of the GNU General
% Public License version 2

phi = cumsum(sqrt(phase_var)*randn(1, size(x, 2)));
for k = 1:size(x, 1);
    x(k, :) = x(k, :).*exp(i*phi);
end;
