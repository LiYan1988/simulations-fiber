function qam_error_vs_snr(modulation)
% Plot the symbol and bit error rate as a function of the SNR.
%
% This function is easy to generalize only for QAM modulation that has
% rectangular Voronoi regions with boundaries parallel with either the
% real or imaginary axes.
%
% INPUTS:
%    modulation : Modulation format (QPSK = 1, 16-QAM = 2)
% OUTPUT:
%    none
%
% Pontus Johannisson, 2010-08-13
% This software is distributed under the terms of the GNU General
% Public License version 2

min_SNR_dB = 0;
max_SNR_dB = 30;
Npoints    = 100; % Number of points to plot for

if (modulation == 1); % QPSK
    s.symb(1, :) = [-1 1];
    s.bits(1, :) = [0 0];
    s.voro(1, :) = [-inf 0 0 inf];

    s.symb(2, :) = [-1 -1];
    s.bits(2, :) = [0 1];
    s.voro(2, :) = [-inf 0 -inf 0];

    s.symb(3, :) = [1 -1];
    s.bits(3, :) = [1 1];
    s.voro(3, :) = [0 inf -inf 0];

    s.symb(4, :) = [1 1];
    s.bits(4, :) = [1 0];
    s.voro(4, :) = [0 inf 0 inf];
elseif (modulation == 2); % 16-QAM
    s.symb(1, :) = [-3 -3];
    s.bits(1, :) = [0 0 0 0];
    s.voro(1, :) = [-inf -2 -inf -2];

    s.symb(2, :) = [-3 -1];
    s.bits(2, :) = [0 0 0 1];
    s.voro(2, :) = [-inf -2 -2 0];

    s.symb(3, :) = [-3 1];
    s.bits(3, :) = [0 0 1 1];
    s.voro(3, :) = [-inf -2 0 2];

    s.symb(4, :) = [-3 3];
    s.bits(4, :) = [0 0 1 0];
    s.voro(4, :) = [-inf -2 2 inf];

    s.symb(5, :) = [-1 -3];
    s.bits(5, :) = [0 1 0 0];
    s.voro(5, :) = [-2 0 -inf -2];

    s.symb(6, :) = [-1 -1];
    s.bits(6, :) = [0 1 0 1];
    s.voro(6, :) = [-2 0 -2 0];

    s.symb(7, :) = [-1 1];
    s.bits(7, :) = [0 1 1 1];
    s.voro(7, :) = [-2 0 0 2];

    s.symb(8, :) = [-1 3];
    s.bits(8, :) = [0 1 1 0];
    s.voro(8, :) = [-2 0 2 inf];

    s.symb(9, :) = [1 -3];
    s.bits(9, :) = [1 1 0 0];
    s.voro(9, :) = [0 2 -inf -2];

    s.symb(10,:) = [1 -1];
    s.bits(10,:) = [1 1 0 1];
    s.voro(10,:) = [0 2 -2 0];

    s.symb(11,:) = [1 1];
    s.bits(11,:) = [1 1 1 1];
    s.voro(11,:) = [0 2 0 2];

    s.symb(12,:) = [1 3];
    s.bits(12,:) = [1 1 1 0];
    s.voro(12,:) = [0 2 2 inf];

    s.symb(13,:) = [3 -3];
    s.bits(13,:) = [1 0 0 0];
    s.voro(13,:) = [2 inf -inf -2];

    s.symb(14,:) = [3 -1];
    s.bits(14,:) = [1 0 0 1];
    s.voro(14,:) = [2 inf -2 0];

    s.symb(15,:) = [3 1];
    s.bits(15,:) = [1 0 1 1];
    s.voro(15,:) = [2 inf 0 2];

    s.symb(16,:) = [3 3];
    s.bits(16,:) = [1 0 1 0];
    s.voro(16,:) = [2 inf 2 inf];
end;
N = size(s.symb, 1); % Number of symbols

% The noise per dimension is given by SNR = mean_power/(2 sigma^2),
% i.e., sigma^2 = mean_power/SNR/2
mean_power = mean(sum(s.symb.^2, 2));
SNR_dB_v   = linspace(min_SNR_dB, max_SNR_dB, Npoints);
SNR_v      = 10.^(SNR_dB_v/10);
sigma_v    = sqrt(mean_power./SNR_v/2);

SER = zeros(size(sigma_v));
BER = zeros(size(sigma_v));
for j = 1:length(sigma_v);
    sigma = sigma_v(j);
    SER_sum = 0;
    BER_sum = 0;
    for k = 1:N;
        for l = 1:N;
            if (k ~= l);
                SER_sum = SER_sum + det_err_prob(s, l, k, sigma)/N;
                BER_sum = BER_sum + det_err_prob(s, l, k, sigma)*bit_err_func(s, l, k)/N;
            end;
        end;
    end;
    SER(j) = SER_sum;
    BER(j) = BER_sum;
end;

figure(1); clf;
semilogy(SNR_dB_v, SER, 'r', SNR_dB_v, BER, 'b');
xlabel('SNR [dB]');
legend('SER', 'BER');
grid on;
zoom on;

figure(2); clf;
semilogx(SER, sigma_v.^2/mean_power, 'r', BER, sigma_v.^2/mean_power, 'b');
ylabel('sigma^2/mean power');
legend('SER', 'BER');
grid on;
zoom on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function P = det_err_prob(s, l, k, sigma)
x1 = s.voro(l, 1) - s.symb(k, 1);
x2 = s.voro(l, 2) - s.symb(k, 1);
y1 = s.voro(l, 3) - s.symb(k, 2);
y2 = s.voro(l, 4) - s.symb(k, 2);
P  = 1/4* ...
     (erf(x1/sqrt(2)/sigma) - erf(x2/sqrt(2)/sigma))* ...
     (erf(y1/sqrt(2)/sigma) - erf(y2/sqrt(2)/sigma));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function E = bit_err_func(s, l, k)
% Compute ratio (incorrect #bits)/(total #bits)

e = sum(s.bits(l, :) ~= s.bits(k, :));
E = e/size(s.bits, 2);
