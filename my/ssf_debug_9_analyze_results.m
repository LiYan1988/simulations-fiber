clc;
clear;
close all;

% Calculate EVM and BER

%%
load debug_9.mat

%%
cidx = 6;
tx = param.signal_tx_complex{cidx};
rx = param.signal_rx_complex{cidx};

objfcn = @(v)v(1)+v(2)*rx - tx;

opts = optimoptions(@lsqnonlin,'Display','off');
x0 = (1+1i)*[1;1]; % arbitrary initial guess
[x_sol,resnorm,residuals,exitflag,output] = lsqnonlin(objfcn,x0,[],[],opts);
x_sol,resnorm,exitflag,output.firstorderopt

rx2 = x_sol(1) + x_sol(2)*rx;
figure; 
hold on;
plot(real(rx2), imag(rx2), '.')
plot(real(tx), imag(tx), 'x')

%% SER
tx_unique = unique(tx);
a = abs(rx2 - tx_unique.');
[u, pidx] = min(a, [], 2);
v = abs(rx2 - tx);

% figure; 
% hold on;
% plot(u)
% plot(v)

ser = sum(u<v)/length(u);

%% EVM
a = abs(rx2-tx).^2;
b = abs(tx).^2;
evm = sqrt(mean(a)/mean(b));