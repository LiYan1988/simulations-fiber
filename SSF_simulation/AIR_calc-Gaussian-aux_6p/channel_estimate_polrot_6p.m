function [parametri,h] = channel_estimate_polrot_6p(x,y,parametri_iniziali,eta,N_particle)  %#codegen
%[parametri,h] = channel_estimate_polrot_6p(x,y,parametri_iniziali,eta,N_particle)
%estimates the parameters of the auxiliary channel model that minimizes the mismatch (in the sense of the KL divergence) with
%respect to the real channel
%
% Parameters are updated by a stochastic gradient algorithm aiming at
% minimizing the auxiliary-channel upper bound to the conditional entropy
%
% The considered auxiliary channel is an AWGN channel with arbitrary
% attenuation, phase noise, and polarization noise:
%
% y_k=h_0*R_k*x_k+n_k
%
% con
% x_k: simboli in ingresso (vettore a 2 componenti)
% y_k: simboli in uscita (vettore a 2 componenti)
% n_k: rumore AWGN (vettore a 2 componenti indipendenti, ciascuna con varianza
%      sg2_n)
% h0:  coefficiente di canale reale costante (attenuazione)
% R_k: matrice complessa 2x2 unitaria tempo variante a 4 parametri reali
%      R_k=R_k(a_k,b_k,c_k,d_k)
% a_k,...,d_k: processi di Wiener indipendenti con incrementi Gaussiani con
%              varianze sg2_a,...,sg2_d
%
% R is a generic unitary 2x2 matrix, which can be expressed as
%               [ cos(d)*exp(ib)   sin(d)*exp(ic)]
%  R = exp(ia)* [                                ]
%               [-sin(d)*exp(-ic)  cos(d)exp(-ib)]
%
%
% Input arguments:
% x: input 2D samples
% y: output 2D samples
% parametri_iniziali: initial values of the six parameters of the auxiliary channel model [h0,sg_n,sg_a,sg_b,sg_c,sg_d]
% eta: adaptation coefficients (for each parameter)
% N_particles: number of particles to represent the distributions.

% Parameters to be adjusted for fine tuning of the algorithm
soglia = 0.3;       %threshold below which particle resampling is performed

% Variable inizialization
N_campioni = length(x);

h0  = parametri_iniziali(1);
sgn = parametri_iniziali(2);
sga = parametri_iniziali(3);
sgb = parametri_iniziali(4);
sgc = parametri_iniziali(5);
sgd = parametri_iniziali(6);

lambda = ones(1,N_campioni);
N_eff = zeros(1,N_campioni);
nr = 0;                 % number of resample operations
parametri = zeros(6,N_campioni);

% Particle inizialization
a = rand(N_particle,1)*2*pi;
b = rand(N_particle,1)*2*pi;
c = rand(N_particle,1)*2*pi;
d = rand(N_particle,1)*2*pi;
pesi = (1/N_particle)*ones(N_particle,1);

% Iterative estimation
for k=1:N_campioni
    
    % Genera nuovo set di particles
    da = randn(N_particle,1); %variazioni Gaussiane con varianza unitaria
    db = randn(N_particle,1); %variazioni Gaussiane con varianza unitaria
    dc = randn(N_particle,1); %variazioni Gaussiane con varianza unitaria
    dd = randn(N_particle,1); %variazioni Gaussiane con varianza unitaria
    
    a = a+sga*da;
    b = b+sgb*db;
    c = c+sgc*dc;
    d = d+sgd*dd;

    % Calcola termini della matrice di rotazione
    expia = exp(1i*a);
    expib = exp(1i*b);
    expic = exp(1i*c);
    cosd  = cos(d);
    sind  = sin(d);
    
    R11 = expia.*expib.*cosd;
    R12 = expia.*expic.*sind;
    R21 = -expia.*conj(expic).*sind;
    R22 = expia.*conj(expib).*cosd;   
    
    % Calcola termini in xy
    yk  = y(:,k);
    xk  = x(:,k);
    xyvec = reshape(conj(yk)*transpose(xk),4,1);
        
    % Calcola la pdf condizionata
    yRx  = [R11 R12 R21 R22]*xyvec;
    d2   = h0^2*sum(abs(xk).^2)+sum(abs(yk).^2)-2*h0*real(yRx);
    sg2n = sgn^2;
    pdf  = (1.0/(pi*sg2n))^2*exp(-(1.0/sg2n)*d2);
    
    % Calcola i nuovi pesi non normalizzati
    pesi = pesi.*pdf;
       
    % Calcola il coefficiente di normalizzazione lambda(k)
    lambda(k)=1./sum(pesi,1);
    
    % Calcola le derivate del termine d2 rispetto ai parametri di attenuazione e rotazione
    yRx_der_b  = [1i*R11 R12 R21 -1i*R22]*xyvec;
    yRx_der_c  = [R11 1i*R12 -1i*R21 R22]*xyvec;
    yRx_der_d  = expia.*([-expib.*sind expic.*cosd -conj(expic).*cosd -conj(expib).*sind]*xyvec);
    d2_der_h0  = 2*h0*sum(abs(xk).^2)-2*real(yRx);
    d2_der_sga = 2*h0*da.*imag(yRx);    
    d2_der_sgb = -2*h0*db.*real(yRx_der_b);
    d2_der_sgc = -2*h0*dc.*real(yRx_der_c);
    d2_der_sgd = -2*h0*dd.*real(yRx_der_d);

    % Calcola le derivate di log(1/lambda(k))=log(sum(pesi)) rispetto ai 6 parametri:
    % d(log(sum(pesi))/dx = 1/sum(pesi)*sum(d(pesi)/dx)=lambda*sum(d(pesi)/dx)
    lmbsusg2   = -lambda(k)/sg2n;
    der_h0     = lmbsusg2*sum(pesi.*d2_der_h0);
    der_sgn    = (2*lmbsusg2/sgn)*sum(pesi.*(2*sg2n-d2));
    der_sga    = lmbsusg2*sum(pesi.*d2_der_sga);
    der_sgb    = lmbsusg2*sum(pesi.*d2_der_sgb);
    der_sgc    = lmbsusg2*sum(pesi.*d2_der_sgc);
    der_sgd    = lmbsusg2*sum(pesi.*d2_der_sgd);
    
    % Normalizza i nuovi pesi
    pesi=lambda(k)*pesi;
    
    % Aggiorna i parametri secondo il gradiente
    h0    = h0+eta(1,k)*der_h0;
    sgn   = sgn+eta(2,k)*der_sgn;
    sga   = sga+eta(3,k)*der_sga;
    sgb   = sgb+eta(4,k)*der_sgb;
    sgc   = sgc+eta(5,k)*der_sgc;
    sgd   = sgd+eta(6,k)*der_sgd;
    parametri(:,k) = [h0;sgn;sga;sgb;sgc;sgd];

    % Verifica la qualità delle particles ed eventualmente le ricampiona
    N_eff(k)=1./sum(pesi.^2);
    if (N_eff(k)<soglia*N_particle)
        nr = nr+1;
        theta = [a b c d];
        [theta,pesi]=particle_resample(theta,pesi);
        a = theta(:,1);
        b = theta(:,2);
        c = theta(:,3);
        d = theta(:,4);
    end
    
end

h = cumsum(log2(lambda))./(1:N_campioni);

end