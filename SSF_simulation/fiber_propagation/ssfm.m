%**************************************************************************
% Propagazione del segnale PDM in un tratto di fibra dispersiva e non
% lineare mediante SSFM (split-step Fourier method) simmetrico.
%**************************************************************************
%
% W1(nfft) - double complex array - Polarization 1
%   Input:  componenti frequenziali del segnale all'ingresso della fibra
%   Output: componenti frequenziali del segnale propagato.
%
% W2(nfft) - double complex array - Polarization 2
%   Input:  componenti frequenziali del segnale all'ingresso della fibra
%   Output: componenti frequenziali del segnale propagato.
%
% nfft - integer
%   Input:  numero di armoniche totali del segnale
%
% F(nfft) - double complex array -
%   Input:  funzione di trasferimento di un passo di fibra
%
% Fh(nfft) - double complex array -
%   Input:  funzione di trasferimento di mezzo passo di fibra
%
% gz(NS) - double precision array -
%    Input: coeff. NL eff. normalizzato in ciascun passo di fibra
%
% NS - integer
%    Input: numero di passi per la propagazione SSFM
%
%**************************************************************************
% Note
% 1) Le componenti frequenziali di W1\2, F e Fh devono essere arrangiate
%    nello stesso ordine (qualsiasi) nei corrispondenti vettori.
% 2) gz(1) contiene il coeff. NL eff. normalizzato calcolato nell'
%    intervallo [0,dz], gz(2) nell'intervallo [dz,2*dz] e cosi' via,
%    avendo diviso la fibra in NS intervalli uguali.
% 3) Il primo e l'ultimo passo richiedono che la propagazione lineare
%    sia valutata solo su mezzo passo (SSFM simmetrico)
%**************************************************************************
%
function Wout = ssfm(Win,nfft,F,Fh,gz,NS,pol)

rn = 1/nfft;

%*** Primo passo (mezzo passo)
Wout(:,1) = Win(:,1).*Fh;                 % Filtraggio lineare
if (pol == 2)
    Wout(:,2) = Win(:,2).*Fh;             % Filtraggio lineare
end

Wout(:,1) = nfft*ifft(Wout(:,1),nfft);                  % Torna nel tempo
if (pol == 2)
    Wout(:,2) = nfft*ifft(Wout(:,2),nfft);              % Torna nel tempo
end

%* Applica le nonlinearita
if (pol == 2)
    fnl = -gz(1)*(Wout(:,1).*conj(Wout(:,1))+Wout(:,2).*conj(Wout(:,2)));
else
    fnl = -gz(1)*(Wout(:,1).*conj(Wout(:,1)));
end
ejgz = rn*exp(1i*wrapToPi(fnl));
Wout(:,1) = Wout(:,1).*ejgz;
if (pol == 2)
    Wout(:,2) = Wout(:,2).*ejgz;
end

Wout(:,1) = fft(Wout(:,1),nfft);           % Torna in frequenza
if (pol == 2)
    Wout(:,2) = fft(Wout(:,2),nfft);       % Torna in frequenza
end

%*** Passi intermedi (regolari)
for i=2:NS
    Wout(:,1) = Wout(:,1).*F;              % Filtraggio lineare
    if (pol == 2)
        Wout(:,2) = Wout(:,2).*F;          % Filtraggio lineare
    end
    
    Wout(:,1) = nfft*ifft(Wout(:,1),nfft);       % Torna nel tempo
    if (pol == 2)
        Wout(:,2) = nfft*ifft(Wout(:,2),nfft);   % Torna nel tempo
    end
    
    %* Applica le non linearita
    if (pol == 2)
        fnl = -gz(i)*(Wout(:,1).*conj(Wout(:,1))+Wout(:,2).*conj(Wout(:,2)));
    else
        fnl = -gz(i)*(Wout(:,1).*conj(Wout(:,1)));
    end
    ejgz = rn*exp(1i*wrapToPi(fnl));
    Wout(:,1) = Wout(:,1).*ejgz;
    if (pol == 2)
        Wout(:,2) = Wout(:,2).*ejgz;
    end
    Wout(:,1) = fft(Wout(:,1),nfft);       % Torna in frequenza
    if (pol == 2)
        Wout(:,2) = fft(Wout(:,2),nfft);   % Torna in frequenza
    end
end

%*** Ultimo passo (mezzo passo, solo dispersione)
Wout(:,1) = Wout(:,1).*Fh;      			 % Filtraggio lineare
if (pol == 2)
    Wout(:,2) = Wout(:,2).*Fh;      		 % Filtraggio lineare
end
