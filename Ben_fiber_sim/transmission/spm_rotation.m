function u = spm_rotation(p, u)
% Function to perform a power-dependent phase rotation modelling the
% distortion caused by self-phase modulation in a fiber.
% The model assumes dispersionless fiber.
% Inputs:
%   p: Parameter structure
%     p.fiber.gamma: Fiber nonlinear coefficient [1/W/m]
%     p.fiber.L    : Fiber length [m]
%     p.fiber.alpha: Fiber attenuation [Np/m]
%   u: Electric field of optical signal
%
% Outputs:
%   u: Electric field of optical signal
%
% Benjamin Foo, Chalmers University of Technology, 2018-08-10

if p.fiber.alpha == 0
    L_eff = p.fiber.L;
else
    L_eff = (1-exp(-p.fiber.alpha*p.fiber.L))/(p.fiber.alpha); % Nonlinear effective length
end

u = u.*exp(1j*p.fiber.gamma*L_eff*abs(u).^2);
end