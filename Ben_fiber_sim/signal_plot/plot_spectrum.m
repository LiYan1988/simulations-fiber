function plot_spectrum(p, u, plot_opts, logscale, esd)
% Plot the spectrum of the input signal.
%
% INPUTS:
%    p         : The parameter struct
%       p.dt_samp : Sampling interval [s]
%       p.f       : Frequency vector [Hz]
%    u         : The signal struct
%    plot_opts : Options to plot function
%    logscale  : Use dB on vertical axis
%    esd       : Plot the energy spectral density (not PSD)
% OUTPUTS:
%   none
%
% By default, fft is used to plot the PSD. The fifth argument can
% be used to obtain the energy spectral density (ESD) instead. The
% difference is a rescaling which affects the dimensions.
%
% Pontus Johannisson, 2009-05-08
% This software is distributed under the terms of the GNU General
% Public License version 2

if nargin < 3
    plot_opts = '';
end

if nargin < 4
    logscale = 1;
end

if nargin < 5
    esd = 0;
end

if esd
    % The energy spectral density
    ESD = abs(p.dt_samp*fft(u, [], 2)).^2*1e9/1e-9; % [nJ/GHz]
    if logscale;
        plot(fftshift(p.f)/1e9, 10*log10(fftshift(ESD, 2)), plot_opts);
        ylabel('ESD [10 log_{10}(nJ/GHz)]');
    else
        plot(fftshift(p.f)/1e9, fftshift(ESD, 2), plot_opts);
        ylabel('ESD [nJ/GHz]');
    end;
    xlabel('frequency [GHz]');
else
    % The power spectral density
    PSD = abs(sqrt(p.dt_samp/length(p.f))*fft(u, [], 2)).^2/1e-3*1e9; % [mW/GHz]
    if logscale;
        plot(fftshift(p.f)/1e9, 10*log10(fftshift(PSD, 2)), plot_opts);
        ylabel('PSD [10 log_{10}(mW/GHz)]');
    else
        plot(fftshift(p.f)/1e9, fftshift(PSD, 2), plot_opts);
        ylabel('PSD [mW/GHz]');
    end;
    xlabel('frequency [GHz]');
end;
grid on;
zoom on;
box on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return
% The following should be the same
if esd
    trapz(p.t, abs(u).^2, 2) % Energy in time domain
    trapz(fftshift(p.f), fftshift(ESD, 2))/(1e9/1e-9) % Energy in frequency domain
else
    trapz(p.t, abs(u).^2, 2)/(p.dt_samp*length(p.f)) % Mean power in time domain
    trapz(fftshift(p.f), fftshift(PSD, 2), 2)/(1/1e-3*1e9) % Mean power in frequency domain
end
