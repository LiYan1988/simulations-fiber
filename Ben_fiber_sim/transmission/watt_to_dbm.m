function P_dBm = watt_to_dbm(P_Watt)
% Convert the power given in Watt to dBm.
%
% INPUTS:
%    P_Watt : The power [Watt]
% OUTPUTS:
%    P_dBm  : The power [dBm]
% This software is distributed under the terms of the GNU General
% Public License version 2

P_dBm = 10*log10(P_Watt/1e-3);

return
% Test case
watt_to_dbm([0.0005 0.001 0.002 0.01])
% Should be (approx) [-3 0 3 10]
