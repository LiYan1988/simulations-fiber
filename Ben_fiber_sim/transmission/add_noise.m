function u = add_noise(u, P_dbm)
% Add white Gaussian noise to a signal
%
% INPUTS:
%    u       : The input signal
%    P_dbm   : Added noise power [dBm]
%
% OUTPUTS:
%    u       : The noise-loaded signal
%
% The sampling interval must be specified since the noise spectral
% density is constant and the amount of noise that should be added
% depends on the bandwidth of the signal.
%
% Benjamin Foo, 2018-04-11
% Based on code from Pontus Johannisson
% This software is distributed under the terms of the GNU General
% Public License version 2

P_n = 10^(P_dbm/10-3);  % Converting noise power from dBm to W

% Add AWGN, the complex noise amplitude should be such that the noise
% power is P_n
u_n = sqrt(P_n)/sqrt(2)*(randn(size(u)) + 1i*randn(size(u)));
%P_n, mean(abs(u_n).^2), pause % The values are similar
u = u + u_n;

