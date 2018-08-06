%***********************************************************************
%**** SUBROUTINE: CPE
%****
%**** Variabili I/O:
%**** - z:       vettore dei campioni equalizzati in uscita
%**** - y:       vettore dei campioni ricevuti.
%**** - x:	     vettore dei simboli trasmessi.
%**** - Nst:     indice a cui iniziare
%**** - N:       numero di campioni su cui operare
%**** - Nsymb:   memoria CPE
%***********************************************************************


function z = cpe(y,x,Nst,N,Nsymb)

gkn = zeros(1,N);

z = y;

for i=Nst:Nst+N-1   
   gk = (x(i-Nsymb:i-1).')*conj(y(i-Nsymb:i-1))+(x(i+1:i+Nsymb).')*conj(y(i+1:i+Nsymb));   
   gkn(i-Nst+1) = gk/abs(gk); 
   z(i) = gkn(i-Nst+1)*y(i);
end
