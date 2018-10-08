function u = nlse_manakov_tod_raman_noise_GPU(p, u)
% Solve the Manakov equation using the split-step Fourier method.
%
% Notation is in accordance with Agrawal's book "Nonlinear fiber
% optics".
%
% INPUTS:
%    p : The parameter struct
%       p.omega         : Angular frequency vector
%       p.steps_per_L_NL: Number of split steps
%       p.fiber.L       : Propagation distance [m]
%       p.fiber.delta_z : Step size [m]
%       p.fiber.beta2   : GVD parameter [s^2/m]
%       p.fiber.beta3   : TOD parameter [s^3/m]
%       p.fiber.gamma   : Kerr nonlinearity [1/W/m]
%       p.fiber.alpha   : Attenuation [1/m]
%       p.raman.P_pump  : Raman backwards pump power [W]
%       p.raman.alpha   : Attenuation at Raman pump wavelength [1/m]
%       p.raman.gain_coeff    : Raman gain coefficient [1/W/m]
%    u : Input amplitude vector (2xN) [sqrt(W)]
% OUTPUT:
%    u : Output amplitude vector (2xN) [sqrt(W)]
%
% All parameters are assumed to be scalar. Third-order dispersion is
% not included.
%
% Pontus Johannisson, 2012-08-28
% Modified for distributed Raman amplification, Henrik Eliasson 2017-05-15
% Modified to run on GPU, Henrik Eliasson 2017-05-17


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
delta_z_GPU = gpuArray(single(delta_z));

% Define operator for the linear part
D     = (1i*p.fiber.beta2/2*(2*pi*p.f).^2 + 1i*p.fiber.beta3/6*(2*pi*p.f).^3 - p.fiber.alpha/2)*delta_z;
exp_D = exp(D);
exp_D_GPU = gpuArray(single(exp_D));

% Define operator for the Raman gain and the Raman noise
R_pump = p.raman.P_pump.*exp(-p.raman.alpha.*(p.fiber.L-linspace(0,p.fiber.L,p.steps_per_L_NL)));
exp_R = exp(1/2*p.raman.gain_coeff.*p.raman.P_pump.*exp(-p.raman.alpha.*(p.fiber.L-linspace(0,p.fiber.L,p.steps_per_L_NL)))*delta_z);
exp_R_GPU = gpuArray(single(exp_R));
R_ASE = p.raman.nsp.*p.const.h.*p.const.nu.*p.raman.gain_coeff.*R_pump.*delta_z;
P_ASE_GPU = gpuArray(single(sqrt(R_ASE*p.f_symb*p.samp_per_symb/2)));

% Store gamma constant in GPU memory
gamma_GPU = gpuArray(single(p.fiber.gamma));

% Initialize power vector on GPU
P_GPU = gpuArray(zeros(1,length(u(1,:)),'single'))

% We don't want to take three steps (two half, one full) in the loop.
% Compare with the last rows of this fcn.
ux_GPU = gpuArray(single(u(1,:)));
uy_GPU = gpuArray(single(u(2,:)));
ux_GPU = ifft(sqrt(exp_D_GPU).*fft(ux_GPU)); % The first half step
uy_GPU = ifft(sqrt(exp_D_GPU).*fft(uy_GPU));

for k = 1:n
    
    % Add noise due to Raman amplification, consider doing this more seldomly
    if p.raman.pump_factor == 1
        ux_GPU =ux_GPU+P_ASE_GPU(k).*randn(1,length(ux_GPU),'single','gpuArray')+1i*P_ASE_GPU(k)/2.*randn(1,length(ux_GPU),'single','gpuArray');
        uy_GPU =uy_GPU+P_ASE_GPU(k).*randn(1,length(uy_GPU),'single','gpuArray')+1i*P_ASE_GPU(k)/2.*randn(1,length(uy_GPU),'single','gpuArray');
    end
    
    % Calculate power and propagate one step
    P_GPU = abs(ux_GPU).^2 + abs(uy_GPU).^2; % Fastest on GPU
    ux_GPU = ifft(exp_R_GPU(k).*exp_D_GPU.*fft(ux_GPU.*exp(delta_z_GPU*1i*gamma_GPU*P_GPU)));
    uy_GPU = ifft(exp_R_GPU(k).*exp_D_GPU.*fft(uy_GPU.*exp(delta_z_GPU*1i*gamma_GPU*P_GPU)));

end
ux_GPU = ifft(1./sqrt(exp_D_GPU).*fft(ux_GPU)); % The last half step backwards
uy_GPU = ifft(1./sqrt(exp_D_GPU).*fft(uy_GPU));

u = gather([ux_GPU;uy_GPU]);