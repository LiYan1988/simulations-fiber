function param = split_step_single_polarization(param)
% Single palorization split step Fourier

% Store dispersive and nonlinear phase factors to speedup code
param.dispersion2 = exp(0.5*1i*param.beta2*param.f.^2*param.dz); 
param.dispersion3 = exp(1i/6*param.beta3*param.f.^3*param.dz);
param.dispersion = param.dispersion2.*param.dispersion3;
param.hhz = 1i*param.gamma*param.dz_eff; 

% --- Main loop
% scheme: 1/2N -> D -> 1/2N; first half step nonlinear
uu = param.data_mod_t_current;
temp = uu.*exp(abs(uu).^2.*param.hhz/2); % note hhz/2
for n=1:param.zn
    % dispersion
    uu = ift(ft(temp, param.df).*param.dispersion, param.df);
    % nonlinearity and loss
    temp = uu.*exp(abs(uu).^2.*param.hhz).*exp(-0.5*param.alpha*param.dz);
end
uu = temp.*exp(-abs(uu).^2.*param.hhz/2); % Final field

% EDFA reamplifies signals
uu = uu*exp(0.5*param.alpha*param.span_length);

if param.ase_exist
    % Power spectral density (PSD) of ASE noise, G=(exp(alpha*L)-1)*h*nu*nsp
    power_noise = (exp(param.alpha*param.span_length)-1)*param.h*param.nu*param.nsp;
    % Convert PSD to signal power in the time domain
    power_noise = power_noise*param.fn*param.df;
    % add noise
    uu = uu + sqrt(0.5*power_noise)*(randn(size(uu, 1), 1)+1i*randn(size(uu, 1), 1));
end

% temp = fftshift(ifft(uu)).* (param.fn*param.dt)/sqrt(2*pi); % Final spectrum
temp = ft(uu, param.df);


% 
param.data_mod_t_current = uu;
param.data_mod_f_current = temp;
% ---
