clc;
clear;
close all;

% Generate WDM channels with different parameters

%% Parameters

% -------------- Primary parameters
param.fmax = 2*pi*400*1e9; % [Hz]
param.fn = 2^22; % number of spectrum points

param.span_length = 100; % [km], span length
param.beta2 = -2.1683e-05; % [ns^2/km], GVD
param.gamma = 1.27; % [(W*km)^-1], nonlinear coefficient of SMF 
param.alpha = log(10)*0.2/10; % [1/km] in linear, 0.2 dB/km
param.zn = 10000; % number of steps per span

% Channel specific parameters, n channels should have n sets of parameters
param.bandwidth_channel = [32*1e9, 10*1e9]; % [Hz]
param.constellation_size = [16, 2]; % 2=BPSK; 4=QPSK; 8=QAM; 16=16QAM; 32=32QAM; 64=64QAM
param.center_frequency_channel = [0*1e9, 50*1e9]; 

% Filter parameters of each channel, assume raised cosine filters
param.roll_off_filter = [0.1, 0.2]; % Roll-off factor of raised cosine filter
param.symbol_in_filter = [10, 6]; % length of impulse response in symbol
param.shape_filter = {'normal', 'sqrt'}; % shape of raised

% -------------- Secondary parameters
param.df = 2*param.fmax/param.fn; % [Hz], frequency resolution
param.f = param.df*[(0:param.fn/2-1) (-param.fn/2:-1)]';

% Time domain parameters
param.tmax = pi/param.df; % [s]
param.dt = 2*param.tmax/param.fn; 
param.t = param.dt*(-param.fn/2:param.fn/2-1)';

param.dz = param.span_length/param.zn; % [km], step size in z

% Channel specific parameters
param.channel_number = length(param.constellation_size);
param.sample_per_symbol = ceil(param.fmax/pi./param.bandwidth_channel);
param.symbol_number = floor(param.fn./param.sample_per_symbol);
param.bit_per_symbol = log2(param.constellation_size);

% Generate signals
rng(294858)
param.data_mod_t = zeros(size(param.t));

for c = 1:param.channel_number
    % Create filter
    filter_tx = rcosdesign(param.roll_off_filter(c), ...
        param.symbol_in_filter(c), param.sample_per_symbol(c), ...
        param.shape_filter{c});
    % Delay of filter
    delay_filter = param.symbol_in_filter(c)/2*param.sample_per_symbol(c);
    
    % Generate bit stream
    data_bit = randi([0, 1], param.symbol_number(c)-2, param.bit_per_symbol(c));
    
    % Convert bits to integers and then to symbols
    data_symbol = bi2de(data_bit);
    data_mod_t_tmp = qammod(data_symbol, param.constellation_size(c));
    data_mod_t_tmp = modnorm(data_mod_t_tmp, 'avpow', 1)*data_mod_t_tmp;
    
    % Pulse shaping
    % First upsample the symbol stream
    data_mod_t_tmp = upfirdn(data_mod_t_tmp, [1], param.sample_per_symbol(c), 1);
    % Pad zeros at the end
    data_mod_t_tmp = [data_mod_t_tmp; zeros(param.sample_per_symbol(c)-1, 1)];
    % Pass through filter
    data_mod_t_tmp = upfirdn(data_mod_t_tmp, filter_tx, 1, 1); % pulse shaping
    % Overhead needs to be removed at the start and end of signal
    length_overhead = length(data_mod_t_tmp)-length(param.t);
    data_mod_t_tmp = data_mod_t_tmp(delay_filter+1:end-length_overhead+delay_filter);
    data_mod_t_tmp = data_mod_t_tmp.*exp(-1i*2*pi*param.center_frequency_channel(c).*param.t);
    param.data_mod_t = param.data_mod_t+data_mod_t_tmp;
end

param.data_mod_f = fftshift(ifft(param.data_mod_t)).*(param.fn*param.dt)/sqrt(2*pi);
param.f_plot = fftshift(param.f)/(2*pi)/1e9;
figure;
plot(param.f_plot, 10*log10(abs(param.data_mod_f).^2))
grid on;