T = 48;                                 % time window (period)
nt = 2^10;                              % number of points
dt = T/nt;                              % timestep (dt)
t = ((1:nt)'-(nt+1)/2)*dt;              % time vector

dz = 5;                            		% total distance per step

u0 = gaussian(t,0,2*sqrt(log(2)));

betap = [0,0,0,1];						% dispersion polynomial
u1 = sspropc(u0,dt,dz,1,0,betap,0);

betap = [0,0,1,1];						% dispersion polynomial
u2 = sspropc(u0,dt,dz,1,0,betap,0);


plot(t,abs(u0).^2,':',...
     t,abs(u1).^2,'-',...
     t,abs(u2).^2,'--');
xlim([-12 12]);
xlabel ('(t-\beta_1z)/T_0');
ylabel ('|u(z,t)|^2/P_0');
