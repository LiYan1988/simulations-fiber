clc
Fs = 200;           % sampling frequency in Hz (you can change it)
Ts = 1/Fs;          % sampling interval
Tot_time = 100;             % total time in seconds (you can change it)
width = 20;                 % width of NRZ pulse in seconds (you can change it)
time = 0:Ts:Tot_time-Ts;    % time on x-axis
modulous = mod(time, width);   
signal = (2*(modulous<width/2))-1; 
plot(time,signal);                 % plots the signal versus time
axis([0 Tot_time -1.5 1.5]);       % adjusts the range of x-axis and y-axis 