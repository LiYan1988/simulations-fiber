function u = nlse_manakov(p, u)
% Solve the Manakov equation using the split-step Fourier method.
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
%    u : Input amplitude vector (2xN) [sqrt(W)]
% OUTPUT:
%    u : Output amplitude vector (2xN) [sqrt(W)]
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
if isinf(p.fiber.delta_z);
    n = 1;
else
    n = ceil(p.fiber.L/p.fiber.delta_z); % Number of iterations
end;
delta_z = p.fiber.L/n; % delta_z * an integer = L

% Define operator for the linear part
%omega = 2*pi/p.T_period*[(0:p.N_samp/2 - 1) (-p.N_samp/2:-1)];
D     = 1/2*(1i*p.fiber.beta2*p.omega.^2 - p.fiber.alpha)*delta_z;
exp_D = exp(D);

% We don't want to take three steps (two half, one full) in the loop.
% Compare with the last rows of this fcn.
u(1, :) = ifft(exp(0.5*D).*fft(u(1, :))); % The first half step
u(2, :) = ifft(exp(0.5*D).*fft(u(2, :)));
for k = 1:n
    %P = abs(u(1, :)).^2 + abs(u(2, :)).^2; % This is less efficient
    P = u(1, :).*conj(u(1, :)) + u(2, :).*conj(u(2, :));
    u(1, :) = ifft(exp_D.*fft(u(1, :).*exp(delta_z*1i*p.fiber.gamma*P)));
    u(2, :) = ifft(exp_D.*fft(u(2, :).*exp(delta_z*1i*p.fiber.gamma*P)));
end
u(1, :) = ifft(exp(-0.5*D).*fft(u(1, :))); % The last half step backwards
u(2, :) = ifft(exp(-0.5*D).*fft(u(2, :)));
