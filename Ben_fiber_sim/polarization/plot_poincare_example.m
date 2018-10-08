function plot_poincare_example()

% Generate all symbols in PM-QPSK
M_pm = exp(1i*pi/2*[0 0 0 0 1 1 1 1 2 2 2 2 3 3 3 3
                   0 1 2 3 0 1 2 3 0 1 2 3 0 1 2 3])*exp(1i*pi/4)/sqrt(2);

figure(1); clf;
plot_poincare();
plot_stokes(jones2stokes(M_pm));
%print -depsc poincare_pmqpsk.eps

% Generate all symbols in PS-QPSK
M_ps = [1+1i -1+1i -1-1i 1-1i 0    0    0   0
        0    0    0   0   1+1i -1+1i -1-1i 1-1i]/sqrt(2);

figure(2); clf;
plot_poincare();
plot_stokes(jones2stokes(M_ps));
%print -depsc poincare_psqpsk.eps
