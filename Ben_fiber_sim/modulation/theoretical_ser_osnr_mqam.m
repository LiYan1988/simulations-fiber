function theoretical_ser_osnr_mqam(M_v)
% A simple function to calculate the theoretical SER as a function of
% OSNR, Eb/N0, and launch power for M-QAM. The nonlinear effects are
% neglected in this function.
%
% INPUTS:
%    M_v : The M vector, QAM constellation size, QPSK = [4];
%
% Pontus Johannisson, 2011-02-21
% This software is distributed under the terms of the GNU General
% Public License version 2

disp('Dual-polarization transmission is assumed')

% Constants
c      = 2.99792458e8;  % Speed of light [m/s]
h      = 6.626076e-34;  % Planck's constant [Js]
lambda = 1.55e-6;       % Wavelength [m]
nu     = c/lambda;      % Frequency [Hz]

% Conversions
att_DB_to_att = log(10)/10/1e3; % * alpha_dB [dB/km]
nm_to_Hz      = c/lambda^2*1e-9;
BW_ref        = 0.1*nm_to_Hz; % 0.1 nm reference bandwidth

% The SMFs
smf.D        = 16.5; % Group-velocity dispersion [ps/(nm km)]
smf.alpha_dB = 0.20; % Attenuation [dB/km]
smf.L        = 80e3; % Transmission distance [m]

% The DCFs
dcf.D        = -120; % Group-velocity dispersion [ps/(nm km)]
dcf.alpha_dB = 0.60; % Attenuation [dB/km]
dcf.L        = -smf.D/dcf.D * smf.L; % Transmission distance [m]

% Amplifier parameters
edfa.G_dB = (smf.alpha_dB*smf.L + dcf.alpha_dB*dcf.L)/1e3; % Amplifier gain [dB]
edfa.G    = 10^(edfa.G_dB/10);                             % Amplifier gain
edfa.n_sp = 1.5;                                           % Spontaneous emission factor

N_spans      = 22;                         % Total number of spans
R_symb       = 14e9;                       % Symbol rate
P_launch_dBm = linspace(-20, 15);          % Launch power per polarization [dBm]
P_launch     = 1e-3*10.^(P_launch_dBm/10); % Launch power per polarization [W]

% PSD of the noise
S_sp = N_spans*edfa.n_sp*h*nu*(edfa.G - 1); % [W/Hz] per polarization

% Calculate the OSNR
OSNR    = (2*P_launch)/(2*S_sp*BW_ref);
OSNR_dB = 10*log10(OSNR);

figure(1); clf;
figure(2); clf;
figure(3); clf;
figure(4); clf;
for M = M_v;
    % Calculate the Eb/N0 (Eb = _average_ energy per bit)
    % P_signal_x = Eb*R_symb*log2(M) % [J/bit] * [symbol/s] * [bit/symbol] = [J/s] = [W]
    % N0 = S_sp
    % We get
    EbN0 = OSNR*(BW_ref/R_symb/log2(M));
    EbN0_dB = 10*log10(EbN0);

    % Calculate Es/N0 (Es = _average_ energy per symbol)
    EsN0 = OSNR*(BW_ref/R_symb);
    EsN0_dB = 10*log10(EsN0);

    % Proakis (4.3-30)
    Q_arg = sqrt(3*log2(M)/(M - 1)*EbN0);
    SER   = 4*(1 - 1/sqrt(M))*0.5*erfc(Q_arg/sqrt(2)).*(1 - (1 - 1/sqrt(M))*0.5*erfc(Q_arg/sqrt(2)));

    figure(1); semilogy(P_launch_dBm, SER); hold on;
    figure(2); semilogy(OSNR_dB,      SER); hold on;
    figure(3); semilogy(EbN0_dB,      SER); hold on;
    figure(4); semilogy(EsN0_dB,      SER); hold on;
    figure(5); plot(P_launch_dBm, OSNR_dB); hold on;
end;
figure(1);
grid on;
axis([-20 15 1e-6 1])
xlabel('Launch power per polarization [dBm]');
ylabel('SER');

figure(2);
grid on;
axis([10 30 1e-6 1])
xlabel('OSNR [dB]');
ylabel('SER');

figure(3);
grid on;
axis([0 30 1e-6 1])
xlabel('Eb/N0 [dB]');
ylabel('SER');

figure(4);
grid on;
axis([0 30 1e-6 1])
xlabel('Es/N0 [dB]');
ylabel('SER');

figure(5);
grid on;
axis([-20 15 10 60])
xlabel('Launch power [dBm]');
ylabel('OSNR [dB]');
