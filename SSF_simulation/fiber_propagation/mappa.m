%
% Mappatura dei bit sui simboli
%
function mp=mappa(i,imod)
      
if (imod == 3)           %BPSK
   mp = 1-2*i;
elseif (imod == 4)       %QPSK
   mp = pskmod(i,4,0,'gray');
elseif (imod == 5)       %16-QAM
   mp = qammod(i,16,0,'gray')./sqrt(10);
elseif (imod == 6)       %64-QAM
   mp = qammod(i,64,0,'gray')./sqrt(42);
end