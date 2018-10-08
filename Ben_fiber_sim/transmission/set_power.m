function u = set_power(u, P)
% Set the optical power to the specified level
%
% INPUTS:
%    u : The input light signal
%    P : The average power [dBm]
% OUTPUTS:
%    u : The output light signal
%
% If there are no output arguments then the power is printed to stdout
% instead.
%
% Pontus Johannisson, 2010-04-20
% This software is distributed under the terms of the GNU General
% Public License version 2

P_ave = mean(abs(u).^2);
P_req = 1e-3 * 10^(P/10);
u = u*sqrt(P_req/P_ave);
