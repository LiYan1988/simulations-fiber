function P_Watt = dbm_to_watt(P_dBm)
% Convert the power given in dBm to Watt.
%
% INPUTS:
%    P_dBm  : The power [dBm]
% OUTPUTS:
%    P_Watt : The power [Watt]
% This software is distributed under the terms of the GNU General
% Public License version 2

P_Watt = 1e-3 * 10.^(P_dBm/10);

return
% Test case
dbm_to_watt([-3 0 3 10])
% Should be (approx) [0.0005 0.001 0.002 0.01]
