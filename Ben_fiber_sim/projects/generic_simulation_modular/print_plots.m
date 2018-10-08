function [power_sort ber_sort evm_sort] = print_plots(results)
Aux_Matrix =zeros(1,length(results));
for ind =1:length(results)
    Aux_Matrix(ind)= results(ind).power;
end
power=Aux_Matrix;
[power_sort,ind_sort]=sort(power);

for ind =1:length(results)
    Aux_Matrix(ind) = results(ind).ber;
end
ber=Aux_Matrix;
ber_sort= ber(ind_sort);
for ind =1:length(results)
    Aux_Matrix(ind) = results(ind).evm_db;
end
evm=Aux_Matrix;
evm_sort=evm(ind_sort);
% for ind =1:length(results) %Solve later.
%     Aux_Matrix(ind) = results(ind).constellation;
% end
% constellation=Aux_Matrix;
% constellation_sort = constellation(ind_sort,:);


% wdm_power=power_sort+10*log(11); %here you need to put the number of channels in the WDM link

figure(1002); plot(power_sort,evm_sort,'-o','LineWidth',2);title('EVM vs Input Power');xlabel('Input Power per Channel [dBm]'); ylabel('EVM [dB]'); grid on;
%figure(1003); semilogy(power_sort,ber_sort,'-o','LineWidth',2);title('BER vs Input Power');xlabel('Input Power per Channel [dBm]'); ylabel('BER]'); grid on;
% for i= 1:length(wdm_power)
%     figure(i); plot(constellation_sort(i,:),'o','LineWidth',2);title(['Constellation for ' num2str(power_sort(i)) 'dBm input power']);xlabel('in-phase'); ylabel('quadrature'); grid on;
% end

%Plot 2 lines
%semilogy(power1,BER1,'-o',power11,BER11,'-o','LineWidth',2);title('BER vs Power per channel');xlabel('Input power per channel[dBm]'); ylabel('BER'); grid on;legend('Single 16QAM Channel','10 OOK + 16QAM');
%plot(power1,EVM1,'-o',power11,EVM11,'-o','LineWidth',2);title('EVM vs Power per channel');xlabel('Input power per channel[dBm]'); ylabel('EVM'); grid on;legend('Single 16QAM Channel','10 OOK + 16QAM');

%You need to import p and u0 to plot the spectrum.
%fig_handle= figure(100); plot_spectrum(p,u0);title(['Modulator output optical waveform']);axis(fig_handle.CurrentAxes, [-400, 400, -60, 25]);

