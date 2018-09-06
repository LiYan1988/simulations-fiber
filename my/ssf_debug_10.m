clc;
close all;
clear;

% How to generate OOK signal without carrier suppression?

%% 
t = -2:0.01:2;
x1 = cos(pi*t);
x2 = 0.5*ones(1, length(x1));

figure; 
hold on;
grid on;
box on;
plot(t, x1)
plot(t, x2);

u = 1/2.*(cos(pi*(x1*0.5+0.5) + 1i*cos(pi*x2)));

figure;
hold on;
plot(t, u)