function P = measure_power(u, info)
% Measure the optical power
%
% INPUTS:
%    u    : The input light signal
%    info : Information string (optional)
% OUTPUTS:
%    P    : The average power [dBm]
%
% If there are no output arguments then the power is printed to stdout
% instead.
%
% Pontus Johannisson, 2010-04-20
% This software is distributed under the terms of the GNU General
% Public License version 2

if nargin < 2
    info = '';
end

P_ave     = mean(abs(u).^2, 2);
P_ave_dBm = 10*log10(P_ave/1e-3);
if nargout < 1
    for k = 1:size(P_ave_dBm, 1)
        if nargin < 2
            fprintf('Mean power, pol. %d: %7.3f dBm = %7.3e W\n', k, P_ave_dBm(k), P_ave(k));
        else
            fprintf('Mean power, pol. %d, %s: %7.3f dBm = %7.3e W\n', k, info, P_ave_dBm(k), P_ave(k));
        end
    end
else
    P = P_ave_dBm;
end
