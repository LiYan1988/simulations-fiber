function param = back_propagation_split_step_single_polarization(param)
% Single palorization split step Fourier

if param.alpha>0
    param.dz_eff_bp = (exp(param.alpha*param.dz)-1)/param.alpha; % [km],
else
    param.dz_eff_bp = param.dz; % [km],
end

% Store dispersive and nonlinear phase factors to speedup code
param.dispersion_bp = exp(-0.5*1i*param.beta2*param.f.^2*param.dz); 
param.hhz_bp = -1i*param.gamma*param.dz_eff_bp; 

% --- Main loop
% scheme: 1/2N -> D -> 1/2N; first half step nonlinear
uu = param.data_mod_t_in_bp;
temp = uu.*exp(abs(uu).^2.*param.hhz_bp/2); % note hhz/2
for n=1:param.zn
    % dispersion
    uu = ift(ft(temp, param.df).*param.dispersion_bp, param.df);
    % nonlinearity and loss
    temp = uu.*exp(abs(uu).^2.*param.hhz_bp).*exp(0.5*param.alpha*param.dz);
end
uu = temp.*exp(-abs(uu).^2.*param.hhz_bp/2); % Final field
% temp = fftshift(ifft(uu)).* (param.fn*param.dt)/sqrt(2*pi); % Final spectrum
temp = ft(uu, param.df);

% 
param.data_mod_t_current_bp = uu;
param.data_mod_f_current_bp = temp;
% ---
