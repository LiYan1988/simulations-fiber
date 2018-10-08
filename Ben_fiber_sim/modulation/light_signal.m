function u = light_signal(p,i)
% Generate the electric field corresponding to a laser. This function
% is simplified right now as no RIN or phase drift is included.
%
% INPUTS:
%    p : The parameter struct
%       p.t_samp          : Time vector [s] (1 x M)
%       p.light.power_dBm : Mean power [dBm] (scalar)
% OUTPUTS:
%    u : The generated light signal (1 x M)
% Modified by Diego Villarani, 2018-09-07
% Pontus Johannisson, 2010-04-19
% This software is distributed under the terms of the GNU General
% Public License version 2

P = 1e-3 * 10^(p.power_array(i)/10); % Power [W]
u = sqrt(P)*ones(size(p.t));
