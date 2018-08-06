%& Generazione della sequenza trasmessa con overlap & save
%& Vengono generati Nu campioni, di cui No sovrapposti ai precedenti

function [ak,Xf]=modoes(ak,Nu,No,imod,sg)

xt = zeros(Nu,1);

%* Overlap dei simboli con i precedenti (No<=Nu)
ak(1:No) = ak(Nu-No+1:Nu);
xt(1:No) = ak(Nu-No+1:Nu);

if (imod == 0) 		  %Continuous wave (CW)
   ak(No+1:Nu) = sg;
   xt(No+1:Nu) = ak(No+1:Nu);
elseif (imod == 1) 	  %Gaussian
   ak(No+1:Nu) = gaussgen(Nu-No,sg*sqrt(0.5));
   xt(No+1:Nu) = ak(No+1:Nu);
elseif (imod == 2) 	  %Constant Envelope Uniform Phase
   ak(No+1:Nu) = sg*exp(2*1i*pi*rand(1,Nu-No));
   xt(No+1:Nu) = ak(No+1:Nu);
elseif (imod == 3) 	  %BPSK
   dataIn = randi([0 1],Nu-No,1);
   dataSym = bi2de(dataIn);   
   ak(No+1:Nu) = sg*mappa(dataSym,imod);
   xt(No+1:Nu) = ak(No+1:Nu);
elseif (imod == 4) 		%QPSK
   dataIn = randi([0 1],Nu-No,2);
   dataSym = bi2de(dataIn);
   ak(No+1:Nu) = sg*mappa(dataSym,imod);
   xt(No+1:Nu) = ak(No+1:Nu);
elseif (imod == 5) 		%16-QAM
   dataSym = randi([0 15],Nu-No,1);
   ak(No+1:Nu) = sg*mappa(dataSym,imod);
   xt(No+1:Nu) = ak(No+1:Nu);
elseif (imod == 6) 		%64-QAM
   dataSym = randi([0 63],Nu-No,1);
   ak(No+1:Nu) = sg*mappa(dataSym,imod);
   xt(No+1:Nu) = ak(No+1:Nu);   
end

Xf = 1/Nu*fft(xt,Nu); %Trasformata 

%NB: fattore 1/Nu serve per rendere la trasformata indipendente dalla banda totale considerata