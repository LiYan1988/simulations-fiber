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
param.snr_per_cloud = cell(1, param.channel_number);
% Total SNR
param.snr_total = cell(1, param.channel_number);
% Points (in Nx2 vector) in the constellation
param.signal_received_constellation = cell(1, param.channel_number);

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
    
    signal = zeros(size(xt_dc, 1), 2);
    signal(:, 1) = real(xt_dc);
    signal(:, 2) = imag(xt_dc);
    
    sumd = zeros(param.sample_per_symbol(cidx), 1);
    for nn=1:param.sample_per_symbol(cidx)
        signal_tmp = downsample(signal, param.sample_per_symbol(cidx), nn-1);
        opts = statset('UseParallel', true);
        [~, ~, sumd_tmp] = kmeans(signal_tmp, ...
            param.constellation_size(cidx), ...
            'Display', 'off', 'maxiter', 1000, ...
            'Replicates', 64, 'Options', opts);
        sumd(nn) = mean(sumd_tmp);
    end
    
    [~, kmean_idx] = min(sumd);
    
    % Sample the signal
    signal = downsample(signal, param.sample_per_symbol(cidx), kmean_idx-1);
    
    % Find center of points
    % idx is the index of centers
    % C is the coordinates of centers
    opts = statset('UseParallel', true);
    [idx, C] = kmeans(signal, param.constellation_size(cidx), ...
        'Display', 'off', 'maxiter', 1000, ...
        'Replicates', 64, 'Options', opts);
    
    % center for each point
    centers = C(idx, :);
    % noise associated to each point
    noise_per_point = sqrt(sum((signal - centers).^2, 2));
    
    % noise per constellation point
    noise_per_cloud = zeros(param.constellation_size(cidx), 1);
    for k=1:param.constellation_size(cidx)
        noise_per_cloud(k) = mean(noise_per_point(idx==k));
    end
    
    %%%%%%%%%%%%%%%%%%%%%
    C_receive_amplitude = sqrt(sum(C.^2, 2));
    C_receive_amplitude = C_receive_amplitude./max(C_receive_amplitude);
    
    % Remove noise in amplitude, this part can be automated by smarter
    % algorithm. Hard code for speed now.
    % This is needed because the received constellation cloud centers are
    % sorted by their amplitudes and angles
    if strcmp(param.channel_type{cidx}, 'ook')
        % OOK
        min_C = min(C_receive_amplitude);
        max_C = max(C_receive_amplitude);
        C_receive_amplitude(C_receive_amplitude==min_C) = 0;
        C_receive_amplitude(C_receive_amplitude==max_C) = 1;
    elseif strcmp(param.channel_type{cidx}, '16qam')
        C_receive_amplitude(C_receive_amplitude>0.8723) = 1;
        C_receive_amplitude((C_receive_amplitude<0.8723)&(C_receive_amplitude>0.5395)) = 0.74;
        C_receive_amplitude(C_receive_amplitude<0.5395) = 0.33;
    end
    
    C_receive_angle = angle(C(:, 1)+1i*C(:, 2))/pi;
    
    [~, C_idx] = sortrows([C_receive_amplitude, C_receive_angle], [1, 2]);
    
    snr_per_cloud = sqrt(sum(C(C_idx, :).^2, 2))./noise_per_cloud(C_idx);
    snr_total = mean(sqrt(sum(centers.^2, 2)))/mean(noise_per_point);
    %%%%%%%%%%%%%%%%%%%%%
    
    %     snr_per_cloud = sqrt(sum(C.^2, 2))./noise_per_cloud;
    %     snr_total = mean(sqrt(sum(centers.^2, 2)))/mean(noise_per_point);
    
    % centers of clouds in constellation diagrams
    param.cloud_centers{cidx} = C;
    % SNR per clound
    param.snr_per_cloud{cidx} = snr_per_cloud;
    % Total SNR
    param.snr_total{cidx} = snr_total;
    % Points (in Nx2 vector) in the constellation
    param.signal_received_constellation{cidx} = signal;
end