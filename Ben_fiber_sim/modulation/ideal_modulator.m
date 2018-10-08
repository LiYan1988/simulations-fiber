function u = ideal_modulator(p, u, e)
% This is an ideal amplitude modulator that shapes the input optical
% exactly according to the electrical signal. Note that this such a
% modulator does not exist.
%
% INPUTS:
%    p : The parameter struct
%    u : The input light signal (1 x M)
%    e : The electric driving signal (1 x M)
% OUTPUTS:
%    u : The modulated light signal (1 x M)
%
% This is an UNREALISTIC modulator that is provided to be able to see
% what happens when the electrical signal is allowed to modulate the
% amplitude of the optical signal. Thus, this modulator will output u
% = u0 * H(e), where the multiplicative factor H = e, where e is
% the input electrical signal.
%
% Pontus Johannisson, 2010-04-20
% This software is distributed under the terms of the GNU General
% Public License version 2

u = u.*e;
