
clc
clear all

addpath('fiber_propagation');
addpath('AIR_calc-Gaussian-aux_3p');

%**************************************************************************
%  Parametri di simulazione
%**************************************************************************
param.lambda = 1.55;      %Lunghezza d'onda [um]
param.Bu = 50;            %Banda segnale osservato = spaziatura WDM [Ghz]

param.pol = 2;            %Polarization division multiplexing (1:no, 2:si)

param.Npunti = 6;         %Number of power values to be considered
param.Pu1dB = -12;          %Potenza segnale osservato in SMF [dBm] (singola polarizzazione)
param.Pinc = 1;            %Potenza segnale osservato in SMF - increment
param.Pu2dB = -16;         %Potenza segnale osservato in DCF [dBm]

param.Nch = 5;             %Numero totale canali WDM
param.Nsc = 1;             %Numero di sottoportanti
param.Nu = 4096/param.Nsc; %Number of symbols(deve essere pari) per ogni subcarrier
param.No = 2048/param.Nsc; %Overlapped symbols(deve essere pari) per ogni subcarrier

param.imod = 1;           %Modulation format (0=CW,1=GAUSS,2=Const.Env.Unif.Phase,3=BPSK,4=QPSK,5=16-QAM, 6=64-QAM)
param.imodw = 1;          %Modulation_format_interferenti_WDM (0=CW,1=GAUSS,2=Const.Env.Unif.Phase,3=BPSK,4=QPSK,5=16-QAM, 6=64-QAM)

param.Nspan = 60;         %Numero di tratte
param.dc = 1;             %Inline dispersion compensation (0=no,1=DCF,2=ideal-per-channel)

%* Noise
param.noise = 3;          %Noise (0=no,1=pre-dec.,2=AWGN in,3=AWGN out,4=inline)
param.ase_source = 0;     %ASE source (0=lumped amp.,1=id. distr.amp.)
param.etasp = 1.6;        %Spont. emission factor

%* Distributed Raman gain(dB/km)
param.ramg1 = 0.0;        %1 Amp.
param.ramg2 = 0.0;        %2 Amp.

%* SMF fiber coefficients
param.alpha1 = 0.2;       %Attenuazione SMF [dB/Km]
param.D1 = 17;            %Dispersione SMF [ps/nm/Km]
param.gm1 = 1.27e-3;      %Coefficiente nonlineare SMF [(W*m)^-1]
param.elle1 = 60;         %Lunghezza SMF [Km]

%* DCF fiber coefficients (or per-channel disperion compensation)
param.alpha2 = 0.57;      %Attenuazione DCF [dB/Km]
param.D2 = -100;          %Dispersione DCF [ps/nm/Km]
param.gm2 = 6.5e-3;       %Coefficiente nonlineare DCF [(W*m)^-1]
param.elle2 = param.elle1*abs(param.D1)/abs(param.D2);  %Lunghezza DCF [Km]

%* SSFM accuracy parameters
param.Nz1f = 1000;        %Numero di passi per span SMF (fw.pr.)
param.Nz1b = 100;         %Numero di passi per span SMF (back. prop.)
param.Nz2f = 1000;        %Numero di passi per span DCF (fw.pr.)
param.Nz2b = 100;         %Numero di passi per span DCF (back. prop.)
param.exbwf = 1.1;        %Fattore di sovracampionamento (banda in eccesso) del segnale WDM (fw.pr.)

%* Backpropagation
param.bp = 1;             %Backpropagation (0=EDC only,1=ideal)
param.detbwb = 1.0;       %Fattore di sovracampionamento del segnale osservato (segnale filtrato per bw.pr.)
param.exbwb = 1.5;        %Fattore di sovracampionamento (banda in eccesso) del segnale osservato (bw.pr.)

%* CPE
param.icpe = 0;           %CPE (0=no,1=on)
param.Ncpe = 20;          %Memoria CPE

%* General
param.seed = 0;              %Seme per generatore casuale (se=0,inizializza con time)
param.Nit_mc = 50*param.Nsc; %Numero di iterazioni per signals propagation
fn_sgs_prop = 'vec/sgs_LA_DM_wdm5_60x60_1sc.mat';

%* AIR calc
N_samples = param.Nit_mc*(param.Nu-param.No);
eta = repmat(linspace(0.00001,0.000001,N_samples),3,1); % step-size coefficient for gradient algorithm (adaptive)
N_particle = 1000;
fn_air_res = 'vec/res_LA_DM_wdm5_60x60_1sc.mat';

%* Blocks activation
act_sgs_prop     = 1;         %Signal propagation
act_sgs_save     = 1;         %Save propagated signals
act_air_res_save = 0;         %Save AIR results
power2sim = 1:param.Npunti;
nsc2sim   = 1:param.Nsc;

%--------------------------------------------------------------------------
%  Generate signals and propagate through optical fiber
%--------------------------------------------------------------------------
if (act_sgs_prop == 1)
    sgs = signals_propagation(param);
    if (act_sgs_save == 1 ) %save propagated signals
        save(fn_sgs_prop,'sgs');
    end
end

%--------------------------------------------------------------------------
%  Load saved signals
%--------------------------------------------------------------------------
if (act_sgs_prop == 0)
    load(fn_sgs_prop);
end

%--------------------------------------------------------------------------
%  3p aux - AIR calc
%--------------------------------------------------------------------------
% Auxiliary channel
for ip_ind=1:length(power2sim)
    ip = power2sim(ip_ind);
    for k_ind=1:length(nsc2sim)
        k = nsc2sim(k_ind);
        sgs(ip).subc(k).snr    = sgs(ip).snr0dB(ip);   %[dB]
        sgs(ip).subc(k).sg2_pn = 0.05;
        sgs(ip).subc(k).h0     = sgs(ip).sg;
        sgs(ip).subc(k).sg2_awgn = abs(sgs(ip).subc(k).h0)^2*10^(-0.1*sgs(ip).subc(k).snr)*2;   %varianza totale (entrambe componenti I e Q)
    end
end

% Adaptation parameters for channel estimate
for ip_ind=1:length(power2sim)
    ip = power2sim(ip_ind);
    for k_ind=1:length(nsc2sim)
        k = nsc2sim(k_ind);
        sgs(ip).subc(k).parametri_iniziali = ...
            [sqrt(sgs(ip).subc(k).sg2_awgn);...
            sqrt(sgs(ip).subc(k).sg2_pn)  ;...
            sgs(ip).subc(k).h0];
    end
end

%* AIR calc
for ip_ind=1:length(power2sim)
    ip = power2sim(ip_ind);
    for k_ind=1:length(nsc2sim)
        k = nsc2sim(k_ind);
        x_tmp = reshape(sgs(ip).subc(k).tx(:,:,1),N_samples,1);
        y_tmp = reshape(sgs(ip).subc(k).rx(:,:,1),N_samples,1);
        res(ip).subc(k) = AIR_calc_aux3p(x_tmp(2:end),y_tmp(2:end),N_particle,eta,sgs(ip).subc(k).parametri_iniziali);
    end
end

%* Calcolo AIR. Aux channel AWGN
AIR_awgn = zeros(param.Nsc,length(res));
for ip_ind=1:length(power2sim)
    ip = power2sim(ip_ind);
    for k_ind=1:length(nsc2sim)
        k = nsc2sim(k_ind);
        AIR_awgn(k,ip) = res(ip).subc(k).AWGN.air;
    end
end

%* Calcolo AIR. Aux channel AWGN + WIENER PN
AIR_aux3p = zeros(param.Nsc,length(res));
for ip_ind=1:length(power2sim)
    ip = power2sim(ip_ind);
    for k_ind=1:length(nsc2sim)
        k = nsc2sim(k_ind);
        AIR_aux3p(k,ip) = res(ip).subc(k).air(end);
    end
end

AIR_awgn_avg = mean(AIR_awgn);
AIR_aux3p_avg = mean(AIR_aux3p);

%--------------------------------------------------------------------------
%  Save AIR results
%--------------------------------------------------------------------------
if (act_air_res_save == 1)
    save(fn_air_res,'res','AIR_awgn_avg','AIR_aux3p_avg','AIR_awgn','AIR_aux3p');
end
