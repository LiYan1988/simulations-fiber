function param = generate_signals(param)
% Generate signals
param.df = 2*param.fmax/param.fn; % [Hz], frequency resolution
% param.f = param.df*[(0:param.fn/2-1) (-param.fn/2:-1)]';
param.f = (-param.fn/2:param.fn/2-1)'*param.df;

% Time domain parameters
param.tmax = pi/param.df; % [s]
param.dt = 2*param.tmax/param.fn;
param.t = param.dt*(-param.fn/2:param.fn/2-1)';

param.dz = param.span_length/param.zn; % [km], step size in z
param.dz_eff = (1-exp(-param.alpha*param.dz))/param.alpha; % [km], effective length of one step, used in the nonlinear component of SSF

% Channel specific parameters
param.channel_number = length(param.constellation_size);
param.sample_per_symbol = ceil(param.fmax/pi./param.bandwidth_channel);
param.symbol_number = floor(param.fn./param.sample_per_symbol);
param.bit_per_symbol = log2(param.constellation_size);

% Generate signals
rng(param.random_seed)
param.data_mod_t_in = zeros(size(param.t));

% Record generated signals
param.filter_tx_channel = cell(1, param.channel_number);
param.data_bit_channel = cell(1, param.channel_number);
param.data_symbol_channel = cell(1, param.channel_number);
param.data_mod_t_channel = cell(1, param.channel_number);

for c = 1:param.channel_number
    % Create filter
    filter_tx = rcosdesign(param.roll_off_filter(c), ...
        param.symbol_in_filter(c), param.sample_per_symbol(c), ...
        param.shape_filter{c});
    param.filter_tx_channel{c} = filter_tx;
    
    % Delay of filter
    delay_filter = param.symbol_in_filter(c)/2*param.sample_per_symbol(c);
    param.delay_filter_channel(c) = delay_filter;
    
    % Generate bit stream
    data_bit = randi([0, 1], param.symbol_number(c)-2, param.bit_per_symbol(c));
    param.data_bit_channel{c} = data_bit;
    
    % Convert bits to integers and then to symbols
    data_symbol = bi2de(data_bit);
    param.data_symbol_channel{c} = data_symbol;
    
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
    
    % Assign power to the signal
    power_normalized = norm(data_mod_t_tmp)^2/param.fn; % normalized power 
    data_mod_t_tmp = data_mod_t_tmp*sqrt(param.power_channel_time(c)/power_normalized);
    param.data_mod_t_channel{c} = data_mod_t_tmp;
    
    param.data_mod_t_in = param.data_mod_t_in+data_mod_t_tmp;
end

% Fourier transform of waveform
% param.data_mod_f_in = fftshift(ifft(param.data_mod_t_in))*sqrt(param.fn*param.dt);
param.data_mod_f_in = ft(param.data_mod_t_in, param.fn, param.df);
% f axis for plotting the spectrum, unit is [GHz]
% param.f_plot = fftshift(param.f)/(2*pi)/1e9;
param.f_plot = param.f/(2*pi)/1e9;

% Time axis for plotting the spectrum, unit is [mus], microsecond
param.t_plot = param.t*1e6; 

% current signal
param.data_mod_t_current = param.data_mod_t_in;
param.data_mod_f_current = param.data_mod_f_in;