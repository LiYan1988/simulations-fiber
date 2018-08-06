T = 32;                                 % time window (period)
nt = 2^10;                              % number of points
dt = T/nt;                              % timestep (dt)
t = ((1:nt)'-(nt+1)/2)*dt;              % time vector

dz = 2;                            		% total distance per step
nz = 2;                                 % total number of steps

betap = [0,0,+1];						% dispersion polynomial

u = zeros(length(t),nz+1);

u(:,1) = solitonpulse(t);	            % u0 = sech(t)

for ii = 1:nz,
  u(:,ii+1) = sspropc(u(:,ii),dt,dz,1,0,betap,0);
end

plot(t,abs(u).^2);
xlim([-8 8]);
xlabel ('(t-\beta_1z)/T_0');
ylabel ('|u(z,t)|^2/P_0');
