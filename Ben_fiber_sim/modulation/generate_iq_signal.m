function iq = generate_iq_signal(p)
% Set up realistic I and Q channels by processing and filtering the
% input data stream. The output is intended to be used to drive the
% modulator.
%
% INPUTS:
%    p : The parameter struct
%       p.t_samp   : Sampling instants [s]
%       p.omega    : Corresponding angular vector [rad/s]
%       p.dt_samp  : Sample duration [s]
%       p.dt_symb  : Symbol duration [s]
%       p.T_period : Time "window" [s]
%       p.N_samp   : Number of samples
%       p.samp_per_symb   : Samples per symbol
%       p.bits_per_symbol : Bits per symbol
%       p.data     : Bit stream
%       p.num_symb : Symbol stream
%       p.modulation.iq_filter_type : Type of filter for IQ signal
%       p.modulation.iq_filter_fwhm : Filter FWHM for IQ signal
%       p.fig.iq : Figure for plotting
% OUTPUT:
%    iq : The I and Q channels (complex)
%
% Pontus Johannisson, 2009-04-29
% This software is distributed under the terms of the GNU General
% Public License version 2

qam = qam_symbols();

% Set up I and Q with one sample per symbol slot
switch p.bits_per_symbol
 case 1 % iq is real
  iq = 2*(p.data - 0.5);
 case 2
  iq = qam.qam4(p.num_symb);
 case 4
  iq = qam.qam16(p.num_symb);
 otherwise
  error('Case not implemented')
end

% Duplicate the ideal signals to the sampling rate
iq = kron(iq, ones(1, p.samp_per_symb));

% Shift the iq signal half a sample duration to get it centered
% exactly in the middle of the symbol slot.
iq = 0.5*(iq + circshift(iq, [0 1]));

% Plot (if requested)
if isfield(p.fig, 'iq')
    figure(p.fig.iq); clf;
    ax(1) = subplot(2, 1, 1); plot(p.t_samp/1e-12, real(iq), 'b+-');hold on;
    ax(2) = subplot(2, 1, 2); plot(p.t_samp/1e-12, imag(iq), 'b+-');hold on;
    linkaxes(ax, 'x');
end

% Filter the I and Q signals to get realistic signals
if ~isfield(p.modulation, 'iq_filter_type') % Default
    % The 3 dB bandwidth is B = 1/p.dt_symb. This means we should set the
    % FWHM bandwidth of the filter to 2/p.dt_symb
    filt.f = p.omega/2/pi; filt.f0 = 0; filt.delta_f = 2/p.dt_symb;
    iq = filter_gaussian(filt, iq);
    fprintf(1, '%s: Using Gaussian IQ filtering\n', mfilename);
else
    switch p.modulation.iq_filter_type
     case 'none'
      ;
     case 'gaussian'
      filt.f = p.omega/2/pi; filt.f0 = 0;
      filt.delta_f = p.modulation.iq_filter_fwhm;
      iq = filter_gaussian(filt, iq);
     otherwise
      error('Unimplemented filter type')
    end
    %iq = ifft(fft(iq).*raised_cosine(f, p.dt_symb, 1));fprintf(1, 'Using raised cosine filtering\n')
end

% Plot (if requested)
if isfield(p.fig, 'iq')
    figure(p.fig.iq);
    subplot(2, 1, 1); plot(p.t_samp/1e-12, real(iq), 'r+-')
    grid on;
    ylabel('I channel')
    title(['I and Q channels from ' strrep(mfilename, '_', '\_')])
    subplot(2, 1, 2); plot(p.t_samp/1e-12, imag(iq), 'r+-')
    grid on;
    ylabel('Q channel')
    xlabel('time [ps]')
end
