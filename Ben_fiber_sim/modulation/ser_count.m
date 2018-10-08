function [ser, se] = ser_count(p, sig, data)
% Count the number of symbol errors
%
% INPUTS:
%    p   : The parameter struct
%       p.modformat : Modulation format
%       p.N_symb    : Number of transmitted symbols
%       p.symnum    : Transmitted symbols (1 x N)
%       p.do_plot   : Plot the constellation
%    sig   : The input signal (1 x N)
%    data  : The transmitted symbols (1 x N)
% OUTPUT:
%    ser : The symbol error rate
%    se  : The number of symbol errors
%
% Pontus Johannisson, 2010-01-03
% This software is distributed under the terms of the GNU General
% Public License version 2

if isfield(p, 'do_plot') && p.do_plot;
    figure(); clf;
    plot(u, '.');
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

det_sym = zeros(size(data)); % Detected symbols

if strcmp(p.modformat, 'qam4');
    det_sym(0 < imag(sig)) = 2;
    ind = (0 < real(sig)); det_sym(ind) = det_sym(ind) + 1;
elseif strcmp(p.modformat, 'qam16');
    %det_sym(imag(sig) < -2) = 0;
    det_sym((-2 < imag(sig)) & (imag(sig) < 0)) = 4;
    det_sym(( 0 < imag(sig)) & (imag(sig) < 2)) = 8;
    det_sym(  2 < imag(sig))                    = 12;

    %ind = (real(sig) < -2); det_sym(ind) = det_sym(ind) + 0;
    ind = ((-2 < real(sig)) & (real(sig) < 0)); det_sym(ind) = det_sym(ind) + 1;
    ind = (( 0 < real(sig)) & (real(sig) < 2)); det_sym(ind) = det_sym(ind) + 2;
    ind =  ( 2 < real(sig));                    det_sym(ind) = det_sym(ind) + 3;
else
    error('Unknown modulation format');
end;

se  = sum(data ~= det_sym); % Symbol errors
ser = se/p.N_symb; % Symbol error rate
