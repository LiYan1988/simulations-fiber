function u = edfa(u, G, Fn, nu, delta_t)
% Simulate an Erbium-doped fiber amplifier by multiplying by the gain
% factor and adding white Gaussian noise.
%
% INPUTS:
%    u       : The input signal
%    G       : Amplifier gain [dB]
%    Fn      : Amplifier noise figure [dB]
%    nu      : The frequency [Hz]
%    delta_t : The sampling interval [s]
% OUTPUTS:
%    u       : The amplified signal
%
% For noise-free amplification, either call with only two input
% parameters or set Fn = -inf [dB].
%
% The sampling interval must be specified since the noise spectral
% density is constant and the amount of noise that should be added
% depends on the bandwidth of the signal.
%
% Pontus Johannisson, 2009-05-06
% This software is distributed under the terms of the GNU General
% Public License version 2

G = 10^(G/10); % Power gain from dB
u = sqrt(G)*u; % Noiseless amplification

if ((nargin > 2) && (Fn > -inf));
    % We should now add white Gaussian noise. The power spectral density
    % is (Agrawal, fiber-optic communication systems, 2nd Ed., (8.1.15))
    %
    % S_sp = (G - 1) n_sp h nu
    %
    % Equation (8.1.19) says
    %
    % Fn = 2 n_sp (G - 1)/G,
    %
    % which gives
    %
    % S_sp = G Fn h nu/2
    %
    % This is a one-sided PSD but for the baseband representation we
    % should add this amount of noise for the entire spectrum.

    h   = 6.626076e-34; % Planck's constant [Js]
    Fn  = 10^(Fn/10); % Noise figure from dB
    S   = G*Fn*h*nu/2; % Power spectral density, [W/Hz]
    BW  = 1/delta_t;   % Bandwidth of signal
    P_n = S*BW;        % Noise power

    % Add AWGN, the complex noise amplitude should be such that the noise
    % power is P_n
    u_n = sqrt(P_n)/sqrt(2)*(randn(size(u)) + 1i*randn(size(u)));
    %P_n, mean(abs(u_n).^2), pause % The values are similar
    u = u + u_n;
end;
