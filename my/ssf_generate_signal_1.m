clc;
clear;
close all;

% --- Frequency parameters
% 400 GHz frequency span, frequency axis is [-fmax:fmax-df]/(2*pi)
% fmax is the maximum frequency on one side, the frequency window is 2*fmax
fmax = 2*pi*400*1e9; % [Hz]
% number of spectrum points
fn = 2^22;
% spectrum resolution
df = 2*fmax/fn; % [Hz]
% spectrum axis
% f is used for calculation, for visualization of spectrum use f/(2*pi)
f = df*[(0:fn/2-1) (-fn/2:-1)]';
% ---

% --- Time parameters
% Time domain
% time window size is 2*pi/df, time axis is [-tmax:tmax-dt]
tmax = pi/df; % [s]
% time resolution
dt = 2*tmax/fn; % [s]
% time axis
t = dt*(-fn/2:fn/2-1)';
% ---

% --- Fiber parameters
% D = 17 [ps/nm/km]
% wavelength = 1.55 [um]
% speed_light = 2.99792458e5 [km/s]
% beta2 = -D*wavelength^2/(speed_light*2*pi) [ns^2/km]
distance = 100; % [km], span length
beta2 = -2.1683e-05; % [ns^2/km], GVD
gamma = 1.27; % [(W*km)^-1], nonlinear coefficient of SMF 
alpha = 0.2; % [dB/km], attenuation in dB
alpha = log(10)*alpha/10; % [1/km] in linear
% ---

% --- Spatial parameters
Ld = dt^2/abs(beta2)*1e18; % [km], dispersion distance
% zn = round(20*distance/Ld/gamma^2); % Number of z steps
zn = 1000;
dz = distance/zn; % [km], step size in z
% ---

% --- Signal parameter
bandwidth_channel = 32*1e9; % [Hz], channel bandwidth
samples_per_symbol = ceil(fmax/pi/bandwidth_channel);
symbol_number = floor(fn/samples_per_symbol);

M = 16; % constellation size, 16QAM
bits_per_symbol = 4;
bit_number = symbol_number*bits_per_symbol;
rng(2378567);
% Pad 0 at the start and end of the symbol stream, so that truncation does
% not impact the signal waveform
data_bit = randi([0, 1], symbol_number-2, bits_per_symbol); 
data_symbol = bi2de(data_bit);
data_mod_t = qammod(data_symbol, M)/sqrt(10);
data_mod_t = [0; data_mod_t; 0];

% raised cosine filter
symbols_in_filter = 10; % length of the filter in symbol
raised_cosine_filter = rcosdesign(0.05, symbols_in_filter, ...
    samples_per_symbol, 'normal');

data_mod_t = upfirdn(data_mod_t, [1], samples_per_symbol, 1); % upsample
data_mod_t = [data_mod_t; zeros(samples_per_symbol-1, 1)]; % pad 0 at the end
data_mod_t = upfirdn(data_mod_t, raised_cosine_filter, 1, 1); % pulse shaping
delay_filter = symbols_in_filter/2*samples_per_symbol; % delay of filter
% remove overhead at the start and end of signal
length_overhead = length(data_mod_t)-length(t);
data_mod_t = data_mod_t(delay_filter+1:end-length_overhead+delay_filter); 

% shift signal to center frequency
% data_mod_t = data_mod_t.*exp(-1i*2*pi*100*1e9*t);

% Fourier transform
data_mod_f = fftshift(ifft(data_mod_t)).*(fn*dt)/sqrt(2*pi);

figure;
subplot(2, 1, 1)
plot(t*1e6, abs(data_mod_t).^2)
xlabel('Time (\mus)')
subplot(2, 1, 2)
plot(fftshift(f)/(2*pi)/1e9, 10*log10(abs(data_mod_f).^2))
xlabel('Frequency (GHz)')

% --- Store dispersive phase shifts to speedup code
dispersion = exp(0.5*1i*beta2*f.^2*dz); % phase factor
hhz = 1i*gamma^2*dz; % nonlinear phase factor
% ---

% --- Main loop
% scheme: 1/2N -> D -> 1/2N; first half step nonlinear
uu = data_mod_t;
temp = uu.*exp(abs(uu).^2.*hhz/2); % note hhz/2
for n=1:zn
    f_temp = ifft(temp).*dispersion;
    uu = fft(f_temp);
    temp = uu.*exp(abs(uu).^2.*hhz);
    if mod(n, 10) == 0
        fprintf('Step %d.\n', n);
    end
end
uu = temp.*exp(-abs(uu).^2.*hhz/2); % Final field
temp = fftshift(ifft(uu)).* (fn*dt)/sqrt(2*pi); % Final spectrum
% ---

% ----Plot output pulse shape and spectrum
figure;
subplot(2,1,1)
h3 = plot (t, abs(uu).^2, '--r', 'linewidth', 2, 'displayname', 'Output'); 
subplot(2,1,2)
h4 = plot(fftshift(f)/(2*pi), 10*log10(abs(temp).^2), '--r', 'linewidth', 2, 'displayname', 'Output'); 

%%
figure;
subplot(2, 1, 1)
h6 = plot(fftshift(f)/(2*pi)/1e9, 10*log10(abs(data_mod_f).^2), '--', 'linewidth', 1, 'displayname', 'Input');
subplot(2, 1, 2)
h5 = plot(fftshift(f)/(2*pi)/1e9, 10*log10(abs(temp).^2), '--', 'linewidth', 1, 'displayname', 'Output'); 
legend([h5, h6])
