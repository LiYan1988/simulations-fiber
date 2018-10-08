function u = nlse_spm(p, u)
% Solve the nonlinear Schrödinger equation using the split-step
% Fourier method.
%
% Notation is in accordance with Agrawal's book "Nonlinear fiber
% optics".
%
% INPUTS:
%    p : The parameter struct
%       p.omega         : Angular frequency vector
%       p.fiber.L       : Propagation distance [m]
%       p.fiber.delta_z : Step size [m]
%       p.fiber.beta2   : GVD parameter [s^2/m]
%       p.fiber.gamma   : Kerr nonlinearity [1/W/m]
%       p.fiber.alpha   : Attenuation [1/m]
%    u : Input amplitude vector [sqrt(W)]
% OUTPUT:
%    u : Output amplitude vector [sqrt(W)]
%
% All parameters are assumed to be scalar. Third-order dispersion is
% not included.
%
% Pontus Johannisson, 2009-05-04
% This software is distributed under the terms of the GNU General
% Public License version 2

% Check transmission length
if (p.fiber.L < 0)
    error('Propagation distance is negative');
elseif (p.fiber.L == 0)
    return;
end

% Check if propagation is linear
if isinf(p.fiber.delta_z)
    n = 1;
else
    n = ceil(p.fiber.L/p.fiber.delta_z); % Number of iterations
end
delta_z = p.fiber.L/n; % delta_z * an integer = L

% Make sure beta3 is defined (set to 0 if not specified)
if ~isfield(p.fiber, 'beta3')
    p.fiber.beta3 = 0;
end

write_log(p, now, sprintf('Starting NLSE split-step solver with %d steps (%.2f m per step %d steps per L_NL)', n, delta_z, p.steps_per_L_NL));

% Define operator for the linear part
D     = 1/2*(1i*p.fiber.beta2*p.omega.^2 + 1i/3*p.fiber.beta3.*p.omega.^3 - p.fiber.alpha)*delta_z; % Linear operator including TOD
exp_D = exp(D);

% omega = 2*pi/p.T_period*[(0:p.N_samp/2 - 1) (-p.N_samp/2:-1)];
% D     = 1/2*(1i*p.fiber.beta2*p.omega.^2 - p.fiber.alpha)*delta_z; % Old calculation of D without beta3 term


% We don't want to take three steps (two half, one full) in the loop.
% Compare with the last rows of this fcn.
u = ifft(exp(0.5*D).*fft(u));  % The first half step
for k = 1:n
    fprintf('Print %s of %s\n',num2str(k),num2str(n));
    u = ifft(exp_D.*fft(u.*exp(delta_z*1i*p.fiber.gamma*abs(u).^2)));
end
u = ifft(exp(-0.5*D).*fft(u)); % The last half step backwards
write_log(p, now, 'NLSE solver completed');
