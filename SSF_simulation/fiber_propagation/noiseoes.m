% Aggiorna vettore di N campioni di rumore per Overlap&Save:
% - I primi No campioni del vettore in uscita sono uguali agli ultimi 
%   No del precedente vettore (in ingresso).
% - Gli altri N-No sono i.i.d. gaussiani complessi a simmetria circolare
% - Fornisce in uscita anche la FFT del vettore
%
% Wf(N) - double complex array
%   Output: gli N campioni di rumore nuovi (in frequenza)
%
%
% wt(N) - double complex array
%   Input:  gli N campioni di rumore precedenti (nel tempo)
%   Output: gli N campioni di rumore nuovi (nel tempo)
%
% N - integer
%   Input: numero di campioni
%
% No - integer
%   Input: numero di campioni di overlap
%
% sgn - double precision
%    Input: deviazione standard  di ciascuna componente I o Q
%           dei campioni di rumore (nel tempo)
%
function [wt,Wf]=noiseoes(wt,N,No,sgn)

Nuo = N-No;

rn=1/N;

% Overlap dei campioni con i precedenti (No<=N)
wt(1:No) = wt(N-No+1:N);

% Generazione nuovi campioni di rumore
wt(No+1:N) = sgn*(randn(1,Nuo)+1i*randn(1,Nuo));
   
Wft = rn*wt;
Wf = fft(Wft,N); %Trasformata

