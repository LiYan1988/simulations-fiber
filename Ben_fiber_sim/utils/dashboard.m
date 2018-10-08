%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dashboard for comparing GMI with and without scaling of the receied     %
% signal power.                                                           %
%                                                                         %
% Written by Tobias Eriksson, tobias.eriksson@nokia.com, 20-jun-2016      %
% Edited 23-aug-2017                                                      %
% Original GMI file obtained from http://fehenberger.de/#sourcecode       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables
close all

%Settings
M = 16; %QAM order
SNRdBvect = 45:-2:-15; %SNR span
Nsymbols = 40000; %Number of symbols.

%Init
hDemodQAM=modem.qamdemod('M', M, 'SymbolOrder', 'Gray');
hModQAM=modem.qammod('M',M,'Symbolorder','gray');
rx.con_map=hModQAM.SymbolMapping;
rx.constellation=hModQAM.Constellation;


for k = 1:length(SNRdBvect)
    SNRdB = SNRdBvect(k);
    display(['Starting SNR = ' num2str(SNRdB) ' dB (' num2str(k) '/' num2str(length(SNRdBvect)) ')']);
    
    %Transmitter:
    N = randi(M,1,Nsymbols);    
    x = modulate(hModQAM,N-1);

    %Channel
    y = awgn(x, SNRdB,'Measured');    
       
    %Calculate GMI with scaling of y+noise:
    gmi_scaled(k) = calcGMI_withNormalization(x,y,'Gray');
    
    %Display results:
%     display(sprintf('SNR: %2.2f dB  -- GMI %2.2f bit/2D-symbol', SNRdB,gmi_scaled(k)));
    display(sprintf('GMI %2.2f bit/2D-symbol', gmi_scaled(k)));
end
%% Plotting
figure(1)

SNR = linspace(min(SNRdBvect),max(SNRdBvect),10000);
SNR = 10.^(SNR/10);
C = log2(1+SNR);
plot(10*log10(SNR),C,'k-')
hold on
plot(SNRdBvect,gmi_scaled,'b--')

ylim([0 max(gmi_scaled)+1])

legend('Capacity','GMI')
xlabel('SNR [dB]')
ylabel('GMI [bit/2D-symbol]')
grid on
hold off


xlabel('SNR [dB]')
ylabel('GMI [bit/2D-symbol]')
grid on
hold off
