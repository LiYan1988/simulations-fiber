function fundamental_soliton()
% Set up and perform propagation of the fundamental soliton.
%
% Pontus Johannisson, 2009-05-04
% This software is distributed under the terms of the GNU General
% Public License version 2

%% Free parameters;
% Amplitude change factor (fundamental soliton has 1)
amp_ratio = 1;

% Steps per nonlinear length
p.steps_per_L_NL = 200;

% Data rate etc.
p.f_symb         = 40e9; % The symbol rate [Baud]
p.samp_per_symb  = 16;  % Samples per symbol
p.N_symb         = 64;   % Number of symbols

% Input signal
p.t_pulse = 10e-12; % Pulse width [s]

% Transmission parameters
p.fiber.L     = 100e3; % Transmission distance [m]
p.fiber.D     = 17.25; % Dispersion parameter [ps/nm/km]
p.fiber.gamma = 1.5;   % Kerr nonlinearity [1/W/m]
p.fiber.alpha = 0;     % Attenuation [1/m]

%% Derived parameters
p.f_samp  = p.f_symb*p.samp_per_symb; % Sampling rate [Hz]
p.N_samp  = p.N_symb*p.samp_per_symb;
p.dt_symb = 1/p.f_symb;               % Symbol interval [s]
p.dt_samp = 1/p.f_samp;               % Sampling interval [s]

p.T_period = p.N_samp*p.dt_samp;
p.t_samp   = (0:p.dt_samp:(p.N_samp - 1)*p.dt_samp) - p.T_period/2;
p.f        =    1/p.T_period*[(0:p.N_samp/2 - 1) (-p.N_samp/2:-1)];
p.omega    = 2*pi/p.T_period*[(0:p.N_samp/2 - 1) (-p.N_samp/2:-1)];

% Dispersion
c             = 2.99792458e8;                % Speed of light [m/s]
lambda        = 1.55e-6;                     % Wavelength [m]
p.fiber.beta2 = -lambda^2/2/pi/c*(p.fiber.D*1e-6); % Group-velocity dispersion [s^2/m]

%% Simulate the fundamental soliton
% Set up the field
L_D             = abs(p.t_pulse^2/p.fiber.beta2);
A               = sqrt(-p.fiber.beta2/p.fiber.gamma/p.t_pulse^2)*amp_ratio;
L_NL            = 1/A^2/p.fiber.gamma;
p.fiber.delta_z = L_NL/p.steps_per_L_NL;
sig.u0          = A*sech(p.t_samp/p.t_pulse);

% Propagate the field
sig.u = nlse_spm(p, sig.u0);

% Plot
figure(1); clf;
plot(p.t_samp/1e-9, abs(sig.u0), 'b-+', ...
     p.t_samp/1e-9, abs(sig.u),  'r-x')
xlabel('time [ps]')
ylabel('power [W]')
title('nlse\_spm')
fprintf(1, 'Initial and final fields for the scalar NLSE\n')
fprintf(1, 'Dispersion length: %6.1f km\n', L_D/1e3)
fprintf(1, 'Nonlinear length:  %6.1f km\n', L_NL/1e3)

% Propagate the field with the Manakov solver, rotate the
% polarization state
u0  = sig.u0;
phi = pi/3; % Angle can be chosen arbitrarily
sig.u0(1, :) = u0*cos(phi);
sig.u0(2, :) = u0*sin(phi);
sig.u = nlse_manakov(p, sig.u0);

% Plot
figure(2); clf;
plot(p.t_samp/1e-9, abs(sig.u0(1, :)), 'b-+',  ...
     p.t_samp/1e-9, abs(sig.u0(2, :)), 'b--+', ...
     p.t_samp/1e-9, abs(sig.u(1, :)),  'r-x',  ...
     p.t_samp/1e-9, abs(sig.u(2, :)),  'r--x')
xlabel('time [ps]')
ylabel('power [W]')
title('nlse\_manakov')
fprintf(1, 'Initial and final fields for the Manakov equation\n')
fprintf(1, 'Dispersion length: %6.1f km\n', L_D/1e3)
fprintf(1, 'Nonlinear length:  %6.1f km\n', L_NL/1e3)

