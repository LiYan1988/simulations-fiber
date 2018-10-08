function [N_se, ser] = ser_count_16qam(p, u, phi)
% Count the number of symbol errors in a 16-QAM received signal
%
% Consider this obsolete and replaced by ser_count.m.
%
% INPUTS:
%    p   : The parameter struct
%       p.N_symb    : Number of transmitted symbols
%       p.num_symb  : Transmitted symbols
%       p.phase_rot : Phase rotation
%       p.do_plot   : Plot the constellation
%    u   : The signal for making decisions
%    phi : Phase for constellation rotation
% OUTPUT:
%    N_se : The number of symbol errors
%    ser  : The symbol error rate
%
% Pontus Johannisson, 2009-04-29
% This software is distributed under the terms of the GNU General
% Public License version 2

if nargin < 3
    phi = 0;
end

% The mean amplitude for a 16-QAM constellation is
mean_amp_16qam = (sqrt(2) + 2*sqrt(10) + sqrt(18))/4;

% Normalize the signal
u_norm = mean_amp_16qam/mean(abs(u));
u = u*u_norm;

% Perform the phase rotation
u = u*exp(1i*phi);

if isfield(p, 'do_plot') && p.do_plot;
    figure(); clf;
    plot(u, '+');
    hold on;
    d = 4;
    plot([-d, d], [ 0,  0], 'r--');
    plot([-d, d], [-2, -2], 'r--');
    plot([-d, d], [ 2,  2], 'r--');
    plot([ 0,  0], [-d, d], 'r--');
    plot([-2, -2], [-d, d], 'r--');
    plot([ 2,  2], [-d, d], 'r--');
    axis equal;
end;

% Using conventions from qam_symbols.m
det_sym = zeros(size(p.num_symb));

%det_sym(imag(u) < -2) = 0;
det_sym((-2 < imag(u)) & (imag(u) < 0)) = 4;
det_sym(( 0 < imag(u)) & (imag(u) < 2)) = 8;
det_sym(  2 < imag(u))                  = 12;

%ind = (real(u) < -2); det_sym(ind) = det_sym(ind) + 0;
ind = ((-2 < real(u)) & (real(u) < 0)); det_sym(ind) = det_sym(ind) + 1;
ind = (( 0 < real(u)) & (real(u) < 2)); det_sym(ind) = det_sym(ind) + 2;
ind =  ( 2 < real(u));                  det_sym(ind) = det_sym(ind) + 3;

N_se = sum(p.num_symb ~= det_sym); % Symbol errors
ser  = N_se/p.N_symb; % Symbol error rate
