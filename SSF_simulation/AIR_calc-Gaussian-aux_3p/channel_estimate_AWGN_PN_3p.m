function [parametri,h] = channel_estimate_AWGN_PN_3p(x,y,parametri_iniziali,eta,N_particle)  %#codegen
%[parametri,h] = channel_estimate_AWGN_PN_3p(x,y,parametri_iniziali,eta,N_particle)
%estimates the parameters of the auxiliary channel model that minimizes the mismatch (in the sense of the KL divergence) with
%respect to the real channel
%
% Parameters are updated by a stochastic gradient algorithm aiming at
% minimizing the auxiliary-channel upper bound to the conditional entropy
%
%The considered auxiliary channel is an AWGN channel with arbitrary
%attenuation and Wiener phase noise:
%
%   y_k=h_0*x_k*exp(i*theta_k)+n_k
%where
%   theta_k=theta_{k-1}+sgt*delta_k
%   n_k: complex AWGN samples with variance sgn*2
%   delta_k: real AWGN samples with unitary variance
%The model is characterized by three parameters:
%   h0: complex channel coefficient
%   sgn: standard deviation of AWGN noise
%   sgt: standard deviation of phase variations between samples
%
% Input arguments:
% x: input samples
% y: output samples
% parametri_iniziali: initial values of the three parameters (sgn,sgt,h0)
% eta: adaptation coefficients (for each parameter)
% N_particles: number of particles to represent the distributions.

% Parameters to be adjusted for fine tuning of the algorithm
soglia = 0.3;       %threshold below which particle resampling is performed

% Variable inizialization
N_campioni = length(x);
sgn = parametri_iniziali(1);
sgt = parametri_iniziali(2);
h0  = parametri_iniziali(3);

lambda = ones(1,N_campioni);
N_eff = zeros(1,N_campioni);
nr = 0;                 % number of resample operations
parametri = complex(zeros(3,N_campioni));

% Particle inizialization
theta = rand(N_particle,1)*2*pi;
pesi = (1/N_particle)*ones(N_particle,1);

% Iterative estimation
for k=1:N_campioni
    
    % Genera nuovo set di fasi
    delta=randn(N_particle,1); %variazioni Gaussiane con varianza unitaria
    theta=theta+sgt*delta;
    
    % Calcola la pdf
    sg2n = sgn^2;
    sg3n = sgn^3;
    ejtheta = exp(1i*theta);
    h0xejtheta = h0*x(k)*ejtheta;
    d2 = abs(h0xejtheta-y(k)).^2;
    pdf = 1/(pi*sg2n)*exp(-d2/sg2n);
    
    % Calcola i nuovi pesi non normalizzati
    pesi = pesi.*pdf;
%     pesi(pesi == 0) = 1e-200;
       
    % Calcola il coefficiente di normalizzazione lambda(k)
    lambda(k)=1./sum(pesi,1);
    
    % Calcola le derivate dell'inverso di lambda rispetto ai parametri
    prodim = imag(h0xejtheta*conj(y(k)));
    der_sgn = 2*lambda(k)/sg3n*sum(pesi.*(d2-sg2n));   
    der_sgt = -2*lambda(k)/sg2n*sum(pesi.*delta.*prodim);
    der_h0 = -2*lambda(k)/sg2n*sum(pesi.*(h0*abs(x(k)).^2-real(y(k)*conj(x(k)*ejtheta))));
    
    % Normalizza i nuovi pesi
    pesi=lambda(k)*pesi;
    
    % Aggiorna i parametri secondo il gradiente
    sgn = sgn+eta(1,k)*der_sgn;
    sgt = sgt+eta(2,k)*der_sgt;
    h0 = h0+eta(3,k)*der_h0;
    parametri(:,k) = [sgn;sgt;h0];    

    % Verifica la qualità delle particles ed eventualmente le ricampiona
%    N_eff(k)=1./sumsqr(pesi);
    N_eff(k)=1./sum(pesi.^2);
    if (N_eff(k)<soglia*N_particle)
        nr = nr+1;
        [theta,pesi]=particle_resample(theta,pesi);
    end
    
end

h = cumsum(log2(lambda))./(1:N_campioni);

end