function u = filter_bessel5(p, u)
% Filter the input signal with a Bessel filter of 5th order
%
% INPUTS:
%    p : The parameter struct
%       p.f       : Frequency vector [Hz]
%       p.f0      : Center frequency [Hz]
%       p.delta_f : FWHM bandwidth [Hz]
%    u : Input signal [sqrt(W)]
% OUTPUT:
%    u : Output signal [sqrt(W)]
%
% Pontus Johannisson, 2009-05-27
% This software is distributed under the terms of the GNU General
% Public License version 2

% Code borrowed from Optilux (licensed under GPL) by P. Serena

x = p.f/(p.delta_f/2); % Frequency normalized to the bandwidth

% Bessel 5th order
Bb = 0.3863;
d0 = 945;
d1 = 945;
d2 = 420;
d3 = 105;
d4 = 15;

om  = 2*pi*x*Bb;
om2 = om  .* om;
om3 = om2 .* om;
om4 = om3 .* om;
om5 = om4 .* om;

pre = d0 - d2*om2 + d4*om4;
pim = d1*om - d3*om3 + om5;
Hf  = abs(d0./(pre + 1i*pim));
disp('*** WARNING in filter_bessel5: Removing phase to avoid time shift')
u   = ifft(fft(u) .* Hf);

%figure; plot(p.f, abs(Hf).^2); grid on; zoom on; pause; close;
