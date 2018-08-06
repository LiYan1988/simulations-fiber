T = 100;
nt = 2^12;
dt = T/nt;
t = ((1:nt)'-(nt+1)/2)*dt;
w = wspace(T,nt);
vs = fftshift(w/(2*pi));
z = 0.1;
nz = 2000;
dz = z/nz;

betap = [0,0,0,1];

u0 = gaussian(t,0,2*sqrt(log(2)));

u = sspropc(u0,dt,dz,nz,0,betap,10^2);
U = fftshift(abs(dt*fft(u)/sqrt(2*pi)).^2);

subplot(121);
plot (t,abs(u).^2);
xlim([-2 4]);
grid on;
xlabel ('(t-\beta_1z)/T_0');
ylabel ('|u(z,t)|^2/P_0');

subplot(122);
plot (vs,U);
xlim([-4 4]);
grid on;
xlabel ('(\nu-\nu_0) T_0');
ylabel ('|U(z,\nu)|^2/P_0');
