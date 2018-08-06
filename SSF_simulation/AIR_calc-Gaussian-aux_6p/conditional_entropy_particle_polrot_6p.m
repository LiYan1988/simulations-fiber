function [h,nr] = conditional_entropy_particle_polrot_6p(x,y,N_particle,parametri)
%[h,nr] = conditional_entropy_particle_polrot_6p(x,y,N_particle,parametri)
%Estimate the conditional entropy of the output given the input using the particle
%method and an auxiliary channel with attenuation, AWGN, and phase/polarization noise:
%
%The auxiliary channel is an AWGN channel with arbitrary
%attenuation, phase noise, and polarization noise:
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
%
% Input arguments:
% x: input 2D samples
% y: output 2D samples
% N_particle: number of particles to represent the distributions.
% parametri: values of the six parameters of the auxiliary channel model (h0,sg_n,sg_a,sg_b,sg_c,sg_d]
%
% Output arguments:
% h: conditional entropy
% nr: number of resampling operations

% Parameters to be adjusted for fine tuning of the algorithm
soglia = 0.3;       %threshold below which particle resampling is performed

% Variable inizialization
N_campioni = length(x);

h0  = parametri(1);
sgn = parametri(2);
sga = parametri(3);
sgb = parametri(4);
sgc = parametri(5);
sgd = parametri(6);

ususg2 = 1.0/sgn^2;
usupisg2q = 1.0/(pi*sgn^2)^2;

lambda = ones(1,N_campioni);
N_eff = zeros(1,N_campioni);
nr = 0;                 % number of resample operations

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

    expia = exp(1i*a);
    expib = exp(1i*b);
    expic = exp(1i*c);
    cosd  = cos(d);
    sind  = sin(d);

    % Calcola il campione attenuato e ruotato (per ciascuna particle)
    xa=h0*x(:,k);
    x1ra = expia.*(cosd.*expib.*xa(1)+sind.*expic.*xa(2));
    x2ra = expia.*(-sind.*conj(expic).*xa(1) + cosd.*conj(expib).*xa(2));
    
    % Calcola la pdf condizionata
    d2 = abs(x1ra-y(1,k)).^2+abs(x2ra-y(2,k)).^2;
    pdf = usupisg2q*exp(-ususg2*d2);

    % Calcola i nuovi pesi non normalizzati
    pesi = pesi.*pdf;
       
    % Calcola il coefficiente di normalizzazione lambda(k)
    lambda(k)=1./sum(pesi,1);
 
    % Normalizza i nuovi pesi
    pesi=lambda(k)*pesi;
    
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