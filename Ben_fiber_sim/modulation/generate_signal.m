function u = generate_signal(p)
% Set up an input signal
%
% INPUTS:
%    p : The parameter struct
%       p.t_samp   : Sampling instants [s]
%       p.dt_symb  : Symbol duration [s]
%       p.T_period : Time "window" [s]
%       p.N_symb   : Number of symbols
%       p.data     : Bit stream
%       p.modulation.format         : Modulation format
%       p.modulation.power_norm     : Unity maximum power (boolean, optional)
%       p.modulation.mean_power_dBm : Mean power [dBm]
%       p.modulation.t_pulse        : Pulse width [s]
%       p.fig.u0 : Figure for plotting
% OUTPUT:
%    u : Amplitude vector [sqrt(W)]
%
% Pontus Johannisson, 2009-04-28
% This software is distributed under the terms of the GNU General
% Public License version 2

iq = generate_iq_signal(p);
u  = zeros(size(p.t_samp));
switch p.modulation.format
 case 'cw'
  % Simple continuous wave
  u = ones(size(p.t_samp));
 case 'rz-ook-gauss'
  % RZ-OOK using Gaussian pulses
  for idx = 1:p.N_symb
      t0 = p.modulation.t_pulse;
      u = u + (p.data(idx) == 1)* ...
          (exp(-1/2/t0^2*(-(idx - 1/2)*p.dt_symb              + p.t_samp).^2) + ...
           exp(-1/2/t0^2*(-(idx - 1/2)*p.dt_symb - p.T_period + p.t_samp).^2) + ...
           exp(-1/2/t0^2*(-(idx - 1/2)*p.dt_symb + p.T_period + p.t_samp).^2));
  end
 case 'dpsk_mzm'
  % DPSK using a Mach-Zehnder modulator
  u = iq;
 case 'dpsk_pm'
  % DPSK using a phase modulator
  u = exp(1i*(pi/2 + real(iq)*pi/2));
 case 'qpsk_phase_mod'
  % QPSK using a phase modulator
  u = exp(1i*(real(iq)*pi/4 + imag(iq)*pi/2));
 case 'qpsk_iq_mod'
  % QPSK using an IQ modulator
  u = iq;
 case 'qpsk_iq_mod_asymm'
  % Asymmetric QPSK with A_high/A_low = 3 using an IQ modulator
  iq = iq*exp(-1i*pi/4);
  u  = real(iq) + 1i/3*imag(iq);
 case '16qam'
  % 16-QAM
  u = iq;
 otherwise
  error('Unknown modulation format')
end

% Adjust the power level
if ~isfield(p.modulation, 'power_norm') || ~p.modulation.power_norm
    sig.P0     = abs(u).^2;
    curr_power = sum(sig.P0)/p.N_samp;
    wish_power = 1e-3*10^(p.modulation.mean_power_dBm/10); % dBm -> W
    if curr_power > 0
        u = u/sqrt(curr_power)*sqrt(wish_power);
        sig.P0 = abs(u).^2;
    end
else % Just normalize amplitude
    u = u/max(abs(u));
    sig.P0 = abs(u).^2;
end

% Plot the power (if requested)
if isfield(p.fig, 'u0')
    figure(p.fig.u0); clf;
    ax(1) = subplot(3, 1, 1); semilogy(p.t_samp, sig.P0)
    xlabel('time [s]'); ylabel('power [W]');grid on;
    ax(2) = subplot(3, 1, 2); plot(p.t_samp, real(u))
    xlabel('time [s]'); ylabel('I channel');grid on;
    ax(3) = subplot(3, 1, 3); plot(p.t_samp, imag(u))
    xlabel('time [s]'); ylabel('Q channel');grid on;
    linkaxes(ax, 'x');
end
