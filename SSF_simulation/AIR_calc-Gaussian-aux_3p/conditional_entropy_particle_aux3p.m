function [h,nr] = conditional_entropy_particle_aux3p(x,y,N_particle,parametri)
%[h,nr] = conditional_entropy_particle_aux3p(x,y,N_particle,parametri)
%Estimate the conditional entropy of the output given the input using the particle
%method and an auxiliary channel with attenuation, AWGN, PN.
%
% Input arguments:
% x: input samples
% y: output samples
% N_particle: number of particles to represent the distributions.
% parametri: values of the three parameters of the auxiliary channel model (sgn,sgt,h0)
%
% Output arguments:
% h: conditional entropy
% nr: number of resampling operations

% Parameters to be adjusted for fine tuning of the algorithm
soglia = 0.3;       %threshold below which particle resampling is performed

% Variable inizialization
N_campioni = length(x);

sgn = parametri(1);
sgt = parametri(2);
h0  = parametri(3);
ususg2 = 1.0/sgn^2;
usupisg2 = 1.0/(pi*sgn^2);

lambda = ones(1,N_campioni);
N_eff = zeros(1,N_campioni);
nr = 0;                 % number of resample operations

% Particle inizialization
theta = rand(N_particle,1)*2*pi;
pesi = (1/N_particle)*ones(N_particle,1);

% Iterative estimation
for k=1:N_campioni
    
    % Genera nuovo set di fasi
    delta=randn(N_particle,1); %variazioni Gaussiane con varianza unitaria
    theta=theta+sgt*delta;
    
    % Calcola la pdf
    ejtheta = exp(1i*theta);
    h0xejtheta = (h0*x(k))*ejtheta;
    d2 = abs(h0xejtheta-y(k)).^2;
    pdf = usupisg2*exp(-ususg2*d2);

    % Calcola i nuovi pesi non normalizzati
    pesi = pesi.*pdf;
       
    % Calcola il coefficiente di normalizzazione lambda(k)
    lambda(k)=1./sum(pesi,1);
 
    % Normalizza i nuovi pesi
    pesi=lambda(k)*pesi;
    
    % Verifica la qualità delle particles ed eventualmente le ricampiona
    %% N_eff(k)=1./sumsqr(pesi);
    N_eff(k)=1./sum((pesi).^2);
    if (N_eff(k)<soglia*N_particle)
        nr = nr+1;
        [theta,pesi]=particle_resample(theta,pesi);
    end
    
end

h = cumsum(log2(lambda))./(1:N_campioni);

end