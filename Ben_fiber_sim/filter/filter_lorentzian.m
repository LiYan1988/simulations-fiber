function u = filter_lorentzian(p, u)
% Filter the input signal with a Lorentzian filter
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

df = p.delta_f/2/sqrt(sqrt(2) - 1);
Hf = 1./(1 + (p.f - p.f0).^2/df^2);
u  = ifft(fft(u) .* Hf);

%figure; plot(p.f, abs(Hf).^2); grid on; zoom on; pause; close;
