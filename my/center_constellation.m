function param = center_constellation(param, dispersion_resudual_length)
% Do the following things for all the channels:
% 1. downconvert to the baseband
% 2. compensate for residual dispersion
% 3. obtain constellation diagram
% 4. find centers of point clouds on the constellation diagram

if nargin==1
    dispersion_resudual_length=0;
end

% centers of clouds in constellation diagrams
param.cloud_centers = cell(1, param.channel_number);
% SNR per clound
% param.snr_per_cloud = cell(1, param.channel_number);
% Total SNR
% param.snr_total = cell(1, param.channel_number);
% Points (in Nx2 vector) in the constellation
param.signal_received_constellation = cell(1, param.channel_number);
% received complex signal/symbols
param.signal_rx_complex = cell(1, param.channel_number);
% original transmitted complex signal/symbols
param.signal_tx_complex = cell(1, param.channel_number);

for cidx=1:param.channel_number
    % Downconversion and dispersion compensation
    [xt_dc, ~, ~] = dispersion_compensation(param.data_mod_t_current, cidx,...
        param.beta2, param.beta3, dispersion_resudual_length, param);
    
    %     % Remove head and tail of the signal
    %     signal = zeros(size(xt_dc(param.delay_filter_channel(cidx)+1:...
    %         end-param.delay_filter_channel(cidx)), 1), 2);
    %
    %     % Convert signal from complex to Nx2 vector
    %     signal(:, 1) = real(xt_dc(param.delay_filter_channel(cidx)+1:...
    %         end-param.delay_filter_channel(cidx)));
    %     signal(:, 2) = imag(xt_dc(param.delay_filter_channel(cidx)+1:...
    %         end-param.delay_filter_channel(cidx)));
    
    sym_rm = param.delay_filter_channel(cidx)/param.sample_per_symbol(cidx);
    param.data_mod_t_dc = xt_dc;
    signal = zeros(size(xt_dc, 1), 2);
    signal(:, 1) = real(xt_dc);
    signal(:, 2) = imag(xt_dc);
    
    sumd = zeros(param.sample_per_symbol(cidx), 1);
    q = zeros(param.sample_per_symbol(cidx), 1);
    for nn=1:param.sample_per_symbol(cidx)
        signal_tmp = downsample(signal, param.sample_per_symbol(cidx), nn-1);
        opts = statset('UseParallel', true);
        % remove the first and last several samples 
        [idx_tmp, c_tmp, sumd_tmp] = ...
            kmeans(signal_tmp((sym_rm+1):(param.symbol_number(cidx)-sym_rm), :), ...
            param.constellation_size(cidx), ...
            'Display', 'off', 'maxiter', 1000, ...
            'Replicates', 64, 'Options', opts);
        sumd(nn) = sum(sumd_tmp);
        c_tmp = c_tmp(:, 1) + 1i*c_tmp(:, 2);
        c_tmp = abs(c_tmp-c_tmp.');
        c_tmp(c_tmp==0) = inf;
        c_tmp = min(c_tmp(:));
        q(nn) = c_tmp/sumd(nn);
    end
    
%     [~, kmean_idx] = min(sumd);
    [~, kmean_idx] = max(q);
    
    % Sample the signal
    signal = downsample(signal, param.sample_per_symbol(cidx), kmean_idx-1);
    signal = signal((sym_rm+1):(param.symbol_number(cidx)-sym_rm), :); % remove the first and last several samples
    [signal, signal_received_out, signal_original_out] = ...
        align_signals(signal, sym_rm, param, cidx);
    
    % Find center of points
    % idx is the index of centers
    % C is the coordinates of centers
    opts = statset('UseParallel', true);
    [~, C] = kmeans(signal, param.constellation_size(cidx), ...
        'Display', 'off', 'maxiter', 1000, ...
        'Replicates', 64, 'Options', opts);
    
    % centers of clouds in constellation diagrams
    param.cloud_centers{cidx} = C;
    param.signal_received_constellation{cidx} = signal;
    % Complex rx signal
    param.signal_rx_complex{cidx} = signal_received_out;
    % Complex tx signal
    param.signal_tx_complex{cidx} = signal_original_out;
end