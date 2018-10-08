function u = iq_modulator(p, u, e)
% Modulate the input electric field in an I/Q modulator.
%
% INPUTS:
%    p : The parameter struct
%    u : The input light signal (1 x M)
%    e : The electric driving signal (2 x M)
% OUTPUTS:
%    u : The modulated light signal (1 x M)
%
% An I/Q modulator will output u = u0 * H(e), where the multiplicative
% factor H = 1/2*(cos(w1) + j*cos(w2)), where the angles w1 = pi *
% V_1/V_pi and w2 = pi * V_2/V_pi.  We assume the input electric
% signal is measured in units of V_pi.
%
% Pontus Johannisson, 2010-04-20
% This software is distributed under the terms of the GNU General
% Public License version 2

u = u * 1/2.*(cos(pi*e(1, :)) + 1i*cos(pi*e(2, :)));
