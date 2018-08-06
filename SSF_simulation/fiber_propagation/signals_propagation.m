
function sgs=signals_propagation(param)


%**************************************************************************
%  Costanti
%**************************************************************************
clight = 2.99792458e5;	 %km/s
hplank = 6.62606957e-34; %J*s
%**************************************************************************

%**************************************************************************
%  Parametri di simulazione
%**************************************************************************
lambda = param.lambda;      %Lunghezza d'onda [um]
Bu = param.Bu;              %Banda segnale osservato = spaziatura WDM [Ghz]

pol = param.pol;            %Polarization division multiplexing (1:no, 2:si)

%***************************** For sweeping power for different simulations
Npunti = param.Npunti;      %Number of power values to be considered
Pu1dB =  param.Pu1dB;       %Potenza segnale osservato in SMF [dBm] (singola polarizzazione)
Pinc = param.Pinc;          %Potenza segnale osservato in SMF - increment in [dBm]
Pu2dB = param.Pu2dB;        %Potenza segnale osservato in DCF [dBm]
%*****************************

Nch = param.Nch;            %Numero totale canali WDM
Nsc = param.Nsc;            %Numero di sottoportanti, number of subcarriers
Nu = param.Nu;              %Number of symbols(deve essere pari) per ogni subcarrier
No = param.No;              %Overlapped symbols(deve essere pari) per ogni subcarrier

imod = param.imod;          %Modulation format (0=CW,1=GAUSS,2=Const.Env.Unif.Phase,3=BPSK,4=QPSK,5=16-QAM, 6=64-QAM)
imodw = param.imodw;        %Modulation_format_interferenti_WDM (0=CW,1=GAUSS,2=Const.Env.Unif.Phase,3=BPSK,4=QPSK,5=16-QAM, 6=64-QAM)

Nspan = param.Nspan;        %Numero di tratte
dc = param.dc;              %Inline dispersion compensation (0=no,1=DCF,2=ideal-per-channel)

%* Noise
noise = param.noise;               %Noise (0=no,1=pre-dec.,2=AWGN in,3=AWGN out,4=inline)
ase_source = param.ase_source;     %ASE source (0=lumped amp.,1=id. distr.amp.)
etasp = param.etasp;               %Spont. emission factor

%* Distributed Raman gain(dB/km)
ramg1 = param.ramg1;        %1 Amp.
ramg2 = param.ramg2;        %2 Amp.

%* SMF fiber coefficients
alpha1 = param.alpha1;      %Attenuazione SMF [dB/Km]
D1 = param.D1;              %Dispersione SMF [ps/nm/Km]
gm1 = param.gm1;            %Coefficiente nonlineare SMF [(W*m)^-1]
elle1 = param.elle1;        %Lunghezza SMF [Km], SMF fiber length

%* DCF fiber coefficients (or per-channel disperion compensation)
alpha2 = param.alpha2;      %Attenuazione DCF [dB/Km]
D2 = param.D2;              %Dispersione DCF [ps/nm/Km]
gm2 = param.gm2;            %Coefficiente nonlineare DCF [(W*m)^-1]
elle2 = param.elle2;        %Lunghezza DCF [Km]

%*** SSFM accuracy parameters ***
Nz1f = param.Nz1f;          %Numero di passi per span SMF (fw.pr.) number of steps per SMF span in forward propagation
Nz1b = param.Nz1b;          %Numero di passi per span SMF (back. prop.)
Nz2f = param.Nz2f;          %Numero di passi per span DCF (fw.pr.)
Nz2b = param.Nz2b;          %Numero di passi per span DCF (back. prop.)
exbwf = param.exbwf;        %Fattore di sovracampionamento (banda in eccesso) del segnale WDM (fw.pr.), Excessampling factor (excess band) of the WDM signal (fw.pr.), oversampling factor

%* Backpropagation
bp = param.bp;             %Backpropagation (0=EDC only,1=ideal)
detbwb = param.detbwb;     %Fattore di sovracampionamento del segnale osservato (segnale filtrato per bw.pr.)
exbwb = param.exbwb;       %Fattore di sovracampionamento (banda in eccesso) del segnale osservato (bw.pr.)

%* CPE
icpe = param.icpe;          %CPE (0=no,1=on)
Ncpe = param.Ncpe;          %Memoria CPE

%* Detection
chdet = 1;                  %Canale WDM osservato

%* MC parameters
seed = param.seed;          %Seme per generatore casuale (se=0,inizializza con time)
Nit_mc = param.Nit_mc;      %Numero di iterazioni

%* General parameters
iv = 0;                     %Mode verbose

%**************************************************************************
%*  Inizializza variabili
%**************************************************************************
Nusc = Nu*Nsc; % total number of frequency points in one channel
sg = 1/sqrt(Nsc);
snr0dB = zeros(Npunti,1);

%* Control random number generation
if (seed > 0)
    rng(seed);
else
    rng('shuffle');
end


%**************************************************************************
%  +++ Ciclo sui diversi valori di potenza +++
%**************************************************************************
for ip=1:Npunti
    %***********************************************************************
    % Definizione dei parametri di simulazione
    %***********************************************************************
    
    %** Parametri di propagazione: potenze, attenuazioni, rotazioni, ...
    
    %* Potenza canale osservato all'ingresso di SMF e DCF
    % channel power observed at input of SMF or DCF
    Pu1=10^((Pu1dB)/10);
    Pu2=10^((Pu2dB)/10);
    
    %* Attenuazione totale
    % SMF attenuation
    A1dB  = alpha1; % [dB/Km]
    at1 = log(10)*A1dB/10; %[1/km], convert alpha from dB to linear
    a1 = at1*elle1; % attenuation of SMF fiber
    % DCF attenuation
    A2dB = alpha2;
    at2 = log(10)*A2dB/10;
    a2 = at2*elle2;
    
    %* Dispersione intera tratta (coeff. (2*pi)^2*beta_2*L/2 in [ns^2])
    % Full dispersion
    % GVD in SMF [ns^/km]
    b1 = -D1*lambda^2/(clight*2*pi);
    % GVD in DCF [ns^/km]
    b2 = -D2*lambda^2/(clight*2*pi); %[ns^2/km]
    % full dispersion in SMF
    bt1 = b1*elle1*2*pi^2;
    % full dispersion in DCF
    bt2 = b2*elle2*2*pi^2;
    
    %* Coefficiente non-lineare intera tratta normalizzato a Pu1
    gmp1 = gm1*Pu1; %[1/km]
    gmp2 = gm2*Pu2;
    gm01 = gmp1*elle1;
    gm02 = gmp2*elle2;
    
    %* Lunghezza spezzoni normalizzata forward propagation
    % normalized length of each step in SMF or DCF
    dz1f = 1/Nz1f;
    dz2f = 1/Nz2f;
    
    %* Lunghezza spezzoni normalizzata backward propagation
    % normalized step length per step in backward propagation
    dz1b = 1/Nz1b;
    dz2b = 1/Nz2b;
    
    %* Dispersione per spezzone (coeff. (2*pi)^2*beta_2*dz/2 in [ns^2])
    % dispersion in each step
    %  forward propagation
    bdz1f = bt1*dz1f;
    bdz2f = bt2*dz2f;
    
    %* Dispersione per spezzone (coeff. (2*pi)^2*beta_2*dz/2 in [ns^2])
    %  backward propagation
    bdz1b = -bt1*dz1b;
    bdz2b = -bt2*dz2b;
    
    %* Non-linearita efficace (nel primo/ultimo spezzone per fw./bw.pr.)
    % nonlinearity per step
    %What does these mean?
    if (abs(a1) > 1e-10) % if total power loss is significant
        ge1f = gm01*(1-exp(-a1*dz1f))/a1;
        ge1b = -gm01*(1-exp(-a1*dz1b))/a1;
    else
        ge1f = gm01*dz1f;
        ge1b = -gm01*dz1b;
    end
    if (abs(a2) > 1e-10)
        ge2f = gm02*(1-exp(-a2*dz2f))/a2;
        ge2b = -gm02*(1-exp(-a2*dz2b))/a2;
    else
        ge2f = gm02*dz2f;
        ge2b = -gm02*dz2b;
    end
    
    %* Nonlinear phase-rotation
    if (abs(a1) > 1.e-10)
        fi1 = gm01*(1-exp(-a1))/a1;
    else
        fi1 = gm01;
    end
    if (abs(a2) > 1.e-10)
        fi2 = gm02*(1-exp(-a2))/a2;
    else
        fi2 = gm02;
    end
    fispm = (fi1+fi2)*Nspan; % phase rotation due to SPM? why?
    fixpm = 2*fispm*(Nch-1); % phase rotation due to XPM? why?
    
    %* Parametri dell'espansione: bande, tempi e numero di punti
    % Expansion parameters: bands, times and number of points
    Bw = Nch*Bu;         %Banda totale WDM (Nch canali di banda Bu), total WDM bandwidth
    Nw = Nusc*Nch;       %Numero totale di componenti nella banda del segnale WDM, total number of frequency points
    Ttot = Nusc/Bu;      %Intervallo temporale di espansione, total time window
    df = 1/Ttot;         %Intervallo frequenziale tra le componenti dell'espansione, frequency resolution
    
    %* Numero totale di componenti a frequenza positiva per l'espansione
    %* in serie dei segnali (la banda totale di espansione e' scelta
    %* maggiore della banda WDM di circa un fattore exbwf>=1, in modo da
    %* tener conto dell'aumento di banda dei segnali durante la
    %* propagazione nonlineare), fw.pr. e bw.pr.
    Nt = 2^ceil(log2(exbwf*Nw));     %Numero di componenti espansione (fw.pr.), number of samples in simulation of forward propagation, consider oversampling rate, and make it 2^N
    Nbu = 2^ceil(log2(detbwb*Nusc));   %Numero di componenti espansione (bw.pr.)
    Nb = 2^ceil(log2(exbwb*Nbu));    %Numero di componenti espansione (bw.pr.)
    
    %* Varie grandezze utili nella costruzione e gestione dei vettori
    %* per la propagazione del segnale
    Nchr = (Nch-1)/2; %Numero di canali dx rispetto CUT, number of channels at the right of CUT
    Nchl = Nchr;      %Numero di canali sx rispetto CUT, number of channels at the left of CUT
    if (chdet <= (Nchl+1)) 
        ich = chdet-1;
    else
        ich = chdet-Nch-1;
    end
    Nt2 = Nt/2;
    Nbu2 = Nbu/2;
    Nb2 = Nb/2;
    Nusc2 = Nusc/2;
    Not = No*Nt/Nusc;     %Campioni di overlap per segnale WDM
    No2 = No/2;         %Meta simboli di overlap (valido per ogni ch.)
    
    %* Vettore frequenze, vector frequency
    fnt = [0:df:(Nt2-1)*df -Nt2*df:df:-df]; %Banda totale
    fnb = [0:df:(Nb2-1)*df -Nb2*df:df:-df]; %Banda backward propagation
    fnu = [0:df:(Nusc2-1)*df -Nusc2*df:df:-df]; %Banda segnale osservato CUT
    
    %* Alcuni parametri caratteristici del sistema, total accumulated
    %dispersion
    if (dc==0)
        DLtot = Nspan*D1*elle1;            %Dispersione totale accumulata
        Toff = Nspan*abs(bt1)*Bw/pi;       %Walk-off temporale agli estremi della banda WDM
    else
        DLtot = Nspan*(D1*elle1+D2*elle2); %Dispersione totale accumulata
        Toff = max(abs(Nspan*(bt1+bt2)-bt2),abs(Nspan*(bt1+bt2)-bt1))*Bw/pi; %Walk-off temporale agli estremi della banda WDM
    end
    
    %* Dispersione totale accumulata
    bttot = -DLtot*pi*lambda^2/clight;
    
    %* Deviazione standard di ciascuna componente di rumore (I o Q)
    %* nel tempo (tenendo conto che la potenza di segnale sulla banda
    %* Bu e' unitaria) a seconda del punto di iniezione e deviazione
    %* standard del rumore di fase (incrementi)
    
    %*  Densita' spettrale potenza (normal.) dell'ASE N0 e SNR per canale
    % spectrum power density of ASE and SNR per channel
    if (ase_source == 0)  	   %lumped amplification
        GsuP = (exp(a1)-1)/Pu1;                        %Factor (G-1)/P
        if (dc > 0) 
            GsuP=GsuP+(exp(a2)-1)/Pu2;
        end
        N0 = Nspan*GsuP*hplank*clight/lambda*etasp*1e21;  %Noise PSD normalized to signal power
        snr0 = 1/(N0*Bu);
        snr0dB(ip) = 10*log10(snr0);
    elseif (ase_source == 1)                             %ideal distr. amplification
        GsuP = ramg1*log(10)/10*elle1/Pu1;                %Factor alpha*L/P
        if (dc > 0)
            GsuP = GsuP+ramg2*log(10)/10*elle2/Pu2;
        end
        N0 = Nspan*GsuP*hplank*clight/lambda*etasp*1e21;        %Noise PSD normalized to signal power
        snr0 = 1/(N0*Bu);
        snr0dB(ip) = 10*log10(snr0);
    end
    
    %* Dev. standard (normal.) ASE ingresso e uscita link (banda Bu*Nt/Nu)
    %* e in-line
    sgnt = sqrt(0.5/snr0*Nt/Nusc);       %ASE ingresso link noise=2,3
    sgnd = sqrt(sgnt^2/Nspan);           %ASE inizio di ogni span noise=4
    sgs(ip).snr0dB = snr0dB;
    sgs(ip).sg = sg;
    sgs(ip).Pu = Pu1dB;
    %**********************************************************************
    
    
    % **********************************************************************
    % Visione dati del sistema e controllo
    % **********************************************************************
    if (iv > 0)
        disp('****************  Verifica dati  ******************');
        fprintf('\n');
        fprintf('Banda segnale osservato (GHz)    = %7.1f \n',Bu);
        fprintf('Banda segnale WDM (GHz)          = %7.1f \n',Bw);
        fprintf('Banda totale di espansione (GHz) = %7.1f \n',Nt*df);
        fprintf('Banda espansione backprop. (GHz) = %7.1f \n',Nb*df);
        fprintf('Tempo totale di espansione (ns)  = %7.1f \n',Ttot);
        fprintf('Max. walk-off su banda WDM (ns)  = %7.1f \n',Toff);
        fprintf('Numero di punti FFT (forwprop.)  = %6.0f \n',Nt);
        fprintf('Numero di punti FFT (backprop.)  = %6.0f \n',Nb);
        fprintf('\n');
        fprintf('Dispers. totale accumul. (ps/nm) = %7.1f \n',DLtot);
        fprintf('Rotazione di fase NL da SPM      = %7.2f \n',fispm);
        fprintf('Rotazione di fase NL da XPM      = %7.2f \n',fixpm);
        fprintf('\n');
        disp('***************************************************');
        fprintf('\n');
    end
    
    if (iv > 0)
        disp('Simulazione in corso...');
        fprintf('\n');
    end
    
    % **********************************************************************
    % Costruzione dei vettori per propagazione in fibra
    % construct of vector for propagation
    % **********************************************************************
    %* Dispersione in fibra per ogni passo, fw.pr.: H(fn,dz)
    F1f = fiber(bdz1f,0,fnt);
    F1hf = fiber(0.5*bdz1f,0,fnt);
    F2f = fiber(bdz2f,0,fnt);
    F2hf = fiber(0.5*bdz2f,0,fnt);
    
    %* Modulo di compensazione ideale per canale (fw.pr.)
    if (dc == 2)
        hdcf = fiber(bt2,0,fnt);
    end
    
    %* Dispersione per bw.pr. (segno,banda,step diversi risp. fw.pr.)
    F1b = fiber(bdz1b,0,fnb);
    F1hb = fiber(0.5*bdz1b,0,fnb);
    F2b = fiber(bdz2b,0,fnb);
    F2hb = fiber(0.5*bdz2b,0,fnb);
    
    %* Modulo di compensazione ideale per canale (bw.pr.)
    if (dc == 2)
        hdcb = fiber(-bt2,0,fnb);
    end
    
    %* Non-linearita in fibra per tutti i passi, fw.pr.: geff(zn)*dz*Pu
    g1f = ge1f*exp(-a1*(0:1:Nz1f-1)*dz1f).';
    g2f = ge2f*exp(-a2*(0:1:Nz2f-1)*dz2f).';
    
    %* Non-linearita per bw.pr. (segno,step,atten diversi risp. fw.pr.)
    g1b = ge1b*exp(-a1*(Nz1b-1:-1:0)*dz1b).';
    g2b = ge2b*exp(-a2*(Nz2b-1:-1:0)*dz2b).';
    
    %* Equalizzatore ideale di dispersione residua
    if (bp == 0)
        hdctot = fiber(-bttot,0,fnu);
    end
    %**********************************************************************
    
    
    % *********************************************************************
    % Ciclo di iterazioni per la simulazione MC
    %**********************************************************************
    
    %* Inizializzazione variabili
    Uif = zeros(Nu,Nch,Nsc,pol);
    uit = zeros(Nu,Nch,Nsc,pol);
    
    if (noise == 2 || noise == 3)
        Wnio = zeros(Nt,pol);          %Rumore AGWN-in/out, noise AWGN
        wnio = zeros(Nt,pol);
    elseif (noise == 4)
        Wnln = zeros(Nt,Nspan,pol);    %Rumore in-line
        wnln = zeros(Nt,Nspan,pol);
    end
    
    Wz = zeros(Nt,pol);
    
    Xrpf = zeros(Nusc,pol);
    xrpt = zeros(Nu,Nsc,pol);
    xcmp = zeros(Nu,Nsc,pol);
    
    disold = 0;
       
    %************************************************************************
    % +++ Ciclo iterazioni MC +++
    %************************************************************************
    for it=1:Nit_mc
        %* Generazione segnali trasmessi WDM
        
        %* Generazione simboli canale osservato PDM
        %* (Overlap & Save)
        if (it == 1)
            for isc=1:Nsc
                %* Polarization 1
                [uit(:,1,isc,1),Uif(:,1,isc,1)] = modoes(uit(:,1,isc,1),Nu,0,imod,sg);
                if (pol == 2)
                    %* Polarization 2
                    [uit(:,1,isc,2),Uif(:,1,isc,2)] = modoes(uit(:,1,isc,2),Nu,0,imod,sg);
                end
            end
        else
            for isc=1:Nsc
                %* Polarization 1
                [uit(:,1,isc,1),Uif(:,1,isc,1)] = modoes(uit(:,1,isc,1),Nu,No,imod,sg);
                if (pol == 2)
                    %* Polarization 2
                    [uit(:,1,isc,2),Uif(:,1,isc,2)] = modoes(uit(:,1,isc,2),Nu,No,imod,sg);
                end
            end
        end
        
        %* Generazione simboli canali interferenti PDM
        %* (Overlap & Save)
        for n=2:Nch
            for isc=1:Nsc
                if (it == 1)
                    %* Polarization 1
                    [uit(:,n,isc,1),Uif(:,n,isc,1)] = modoes(uit(:,n,isc,1),Nu,0,imodw,sg);
                    if (pol == 2)
                        %* Polarization 2
                        [uit(:,n,isc,2),Uif(:,n,isc,2)] = modoes(uit(:,n,isc,2),Nu,0,imodw,sg);
                    end
                else
                    %* Polarization 1
                    [uit(:,n,isc,1),Uif(:,n,isc,1)] = modoes(uit(:,n,isc,1),Nu,No,imodw,sg);
                    if (pol == 2)
                        %* Polarization 2
                        [uit(:,n,isc,2),Uif(:,n,isc,2)] = modoes(uit(:,n,isc,2),Nu,No,imodw,sg);
                    end
                end
            end
        end
        
        %uit(symbols , WDM channels , subcariere , Polarizaton)
        %* Generazione segnale multi-carrier
        Uif_mcm = muxmcm(Uif,Nsc,Nu,Nch);%Uf_mux(freqsamples,..., WDM channel, Pol)
        
        %** MUX WDM per la generazione del segnale da propagare Wz,
        %*  esteso con zeri alla banda di fw.pr. Bt
        %* Polarization 1
        Wz(:,1) = muxwdm(Nt,Uif_mcm(:,:,1),Nusc,Nchr,Nchl);
        if (pol == 2)
            %* Polarization 2
            Wz(:,2) = muxwdm(Nt,Uif_mcm(:,:,2),Nusc,Nchr,Nchl);
        end
        
        %* Somma rumore AWGN in ingresso con O&S
        if (noise == 2)
            if (it == 1)
                % Polarization 1
                [wnio(:,1),Wnio(:,1)] = noiseoes(wnio(:,1),Nt,0,sgnt);
                if (pol == 2)
                    % Polarization 2
                    [wnio(:,2),Wnio(:,2)] = noiseoes(wnio(:,2),Nt,0,sgnt);
                end
            else
                % Polarization 1
                [wnio(:,1),Wnio(:,1)] = noiseoes(wnio(:,1),Nt,Not,sgnt);
                if (pol == 2)
                    % Polarization 2
                    [wnio(:,2),Wnio(:,2)] = noiseoes(wnio(:,2),Nt,Not,sgnt);
                end
            end
            % Polarization 1
            Wz(1:Nt,1)=Wz(1:Nt,1)+Wnio(1:Nt,1);
            if (pol == 2)
                % Polarization 2
                Wz(1:Nt,2)=Wz(1:Nt,2)+Wnio(1:Nt,2);
            end
        end
        
        %***   Propagazione in fibra
        for is=1:Nspan
            %* Somma rumore AWGN in-line con O&S
            if (noise == 4)
                if (it == 1)
                    % Polarization 1
                    [wnln(:,is,1),Wnln(:,is,1)] = noiseoes(wnln(:,is,1),Nt,0,sgnd);
                    if (pol == 2)
                        % Polarization 2
                        [wnln(:,is,2),Wnln(:,is,2)] = noiseoes(wnln(:,is,2),Nt,0,sgnd);
                    end
                else
                    % Polarization 1
                    [wnln(:,is,1),Wnln(:,is,1)] = noiseoes(wnln(:,is,1),Nt,Not,sgnd);
                    if (pol == 2)
                        %Polarization 2
                        [wnln(:,is,2),Wnln(:,is,2)] = noiseoes(wnln(:,is,2),Nt,Not,sgnd);
                    end
                end
                % Polarization 1
                Wz(1:Nt,1) = Wz(1:Nt,1) + Wnln(1:Nt,is,1);
                if (pol == 2)
                    % Polarization 2
                    Wz(1:Nt,2) = Wz(1:Nt,2) + Wnln(1:Nt,is,2);
                end
            end
            
            Wz = ssfm(Wz,Nt,F1f,F1hf,g1f,Nz1f,pol);
            
            if (dc == 1)         %fibra DCF
                Wz = ssfm(Wz,Nt,F2f,F2hf,g2f,Nz2f,pol);
            elseif (dc == 2)     %modulo di compensazione per canale
                % Polarization 1
                Wz(:,1) = Wz(:,1).*hdcf;
                if (pol == 2)
                    % Polarization 2
                    Wz(:,2) = Wz(:,2).*hdcf;
                end
            end
        end
        
        %* Somma rumore AWGN in uscita con O&S
        if (noise == 3)
            if (it == 1)
                % Polarization 1
                [wnio(:,1),Wnio(:,1)] = noiseoes(wnio(:,1),Nt,0,sgnt);
                if (pol == 2)
                    % Polarization 2
                    [wnio(:,2),Wnio(:,2)] = noiseoes(wnio(:,2),Nt,0,sgnt);
                end
            else
                % Polarization 1
                [wnio(:,1),Wnio(:,1)] = noiseoes(wnio(:,1),Nt,Not,sgnt);
                if (pol == 2)
                    % Polarization 2
                    [wnio(:,2),Wnio(:,2)] = noiseoes(wnio(:,2),Nt,Not,sgnt);
                end
            end
            % Polarization 1
            Wz(1:Nt,1)=Wz(1:Nt,1)+Wnio(1:Nt,1);
            if (pol == 2)
                % Polarization 2
                Wz(1:Nt,2)=Wz(1:Nt,2)+Wnio(1:Nt,2);
            end
        end
        
        %* Backpropagation: estrae segnale su banda > Bu (Nb punti), aggiunge
        %* rumore ASE e di fase opzionale con O&S e lo
        %* processa su banda Bb>Bu (Nb punti) con SSMF a coeff. opposti
        
        if (bp == 1)
            Ub = zeros(Nb,pol);
            Ub(1:Nbu2,:) = Wz(1:Nbu2,:);
            Ub(Nb-Nbu2+1:Nb,:) = Wz(Nt-Nbu2+1:Nt,:);
            
            for is=1:Nspan
                if (dc == 1)           %fibra DCF
                    Ub = ssfm(Ub,Nb,F2b,F2hb,g2b,Nz2b,pol);
                elseif (dc == 2)       %modulo di compensazione per canale
                    Ub(:,1) = Ub(:,1).*hdcb;
                    if (dual == 2)
                        Ub(:,2) = Ub(:,2).*hdcb;
                    end
                end
                Ub = ssfm(Ub,Nb,F1b,F1hb,g1b,Nz1b,pol);
            end
            
            %* Calcola segnale di uscita uout(t)=ub(t) in
            %* frequenza sulla banda Bu
            Xrpf(1:Nusc2,:) = Ub(1:Nusc2,:);
            Xrpf(Nusc2+1:Nusc,:) = Ub(Nb-Nusc2+1:Nb,:);
            
            %* EDC: Compensazione ideale della dispersione segnale su banda BU
        elseif (bp == 0)
            if (ich > 0)
                Xrpf(1:Nusc2,:) = Wz(ich*Nusc+1:ich*Nusc+Nusc2,:);
                Xrpf(Nusc2+1:Nusc,:) = Wz((ich-1)*Nusc+Nusc2+1:ich*Nusc,:);
            elseif (ich == 0)
                Xrpf(1:Nusc2,:) = Wz(1:Nusc2,:);
                Xrpf(Nusc2+1:Nusc,:) = Wz(Nt-Nusc2+1:Nt,:);
            else
                Xrpf(1:Nusc2,:) = Wz(Nt+ich*Nusc+1:Nt+ich*Nusc+Nusc2,:);
                Xrpf(Nusc2+1:Nusc,:) = Wz(Nt+ich*Nusc-Nusc2+1:Nt+ich*Nusc,:);
            end
            
            Xrpf(:,1) = Xrpf(:,1).*hdctot.*exp(1i*4*pi^2*b1*(fnu.')*Bu*Nspan*elle1*ich);
            if (pol == 2)
                Xrpf(:,2) = Xrpf(:,2).*hdctot.*exp(1i*4*pi^2*b1*(fnu.')*Bu*Nspan*elle1*ich);
            end
        end
        
        Xrpf_demux = demuxmcm(Xrpf,Nsc,Nu,pol);
        
        %* Multicarrier signal ricevuto nel dominio del tempo
        for k=1:pol
            xrpt(:,:,k) = Nu*ifft(squeeze(Xrpf_demux(:,:,k)),Nu).*exp(1i*2*pi^2*b1*Bu^2*Nspan*elle1*abs(ich)^2);
        end
        
        %* Overlap & Save: scarta No/2 campioni iniziali e finale
        diff = abs(xrpt(No2,1)-disold).^2;
        disold = xrpt(Nu-No2,1);
        if ((it > 1) && (iv >0))
            fprintf('Errore stimato Ov&Sav = %12.5e \n \n',diff);
        end
        
        %*** Compensazione opzionale fase
        for k=1:Nsc
            if (icpe == 1)
                for i=1:pol
                   xcmp(:,k,i) = cpe(xrpt(:,k,i),uit(:,chdet,k,i),Ncpe+1,Nu-2*Ncpe,Ncpe);
                end                             
            else
                %*no compensazione
                xcmp = xrpt;
            end
        end
        
        for k=1:Nsc
            sgs(ip).subc(k).tx(:,it,:) = 1/sg*uit(No2+1:Nu-No2,chdet,k,:);
            sgs(ip).subc(k).rx(:,it,:) = xcmp(No2+1:Nu-No2,k,:);
        end
    end
          
    %* Aggiorna parametri per simulazione successiva
    Pu1dB = Pu1dB+Pinc;
    Pu2dB = Pu2dB+Pinc;
    
    disp('done')
end