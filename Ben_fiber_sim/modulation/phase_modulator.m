function u = phase_modulator(p, u, e)
% Modulate the input electric field in a phase modulator.
%
% INPUTS:
%    p : The parameter struct
%    u : The input light signal (1 x M)
%    e : The electric driving signal (1 x M)
% OUTPUTS:
%    u : The modulated light signal (1 x M)
%
% A phase modulator will output u = u0 * H(e), where the
% multiplicative factor H = exp(i * w1), where the angle w1 = pi *
% V_1/V_pi.  We assume the input electric signal is measured in units
% of V_pi.
%
% Pontus Johannisson, 2010-04-20
% This software is distributed under the terms of the GNU General
% Public License version 2

u = u.*exp(1i*pi*e);
