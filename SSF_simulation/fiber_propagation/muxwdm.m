%
% Wz(Nt) - double complex array -
%   In uscita: segnale WDM multiplexato    
%
% Nt - integer
%   In ingresso: numero di armoniche totali per Wz
%
% Uif(Nu,Nch) - double complex array -
%   In ingresso: segnale trasmesso per ogni canale WDM
%
% Nu - integer
%   In ingresso: numero di armoniche totali per ogni canale
%
% Nchr - integer
%   In ingresso: numero di canali a destra del canale osservato
%
% Nchl - integer
%   In ingresso: numero di canali a sinistra del canale osservato
%
% Il segnale WDM viene generato in frequenza rispettando la notazione FFT
% Esempio 5 canali WDM 
% Allocazione canali modello FFT [1,2,3,4,5]
% Rappresentazione segnale 4,5,1,2,3 con 1 canale centrale osservato

function Wz=muxwdm(Nt,Uif,Nu,Nchr,Nchl)

Wz = zeros(Nt,1);
Nu2 = Nu/2;

%* Canale centrale
Wz(1:Nu2) = Uif(1:Nu2,1);
Wz(Nt-Nu2+1:end) = Uif(Nu2+1:end,1);

%* Canali a destra
for i=1:Nchr
    Wz((i-1)*Nu+Nu2+1:i*Nu) = Uif(Nu2+1:Nu,i+1);        %frequenze negative
    Wz(i*Nu+1:i*Nu+Nu2) = Uif(1:Nu2,i+1);               %frequenze positive
end

%* Canali a sinistra
st = Nt-Nu2-Nchl*Nu;
for i=1:Nchl
    Wz(st+(i-1)*Nu+1:st+(i-1)*Nu+Nu2) = Uif(Nu2+1:Nu,Nchr+i+1);   %frequenze negative
    Wz(st+(i-1)*Nu+Nu2+1:st+i*Nu) = Uif(1:Nu2,Nchr+i+1);          %frequenze positive
end





