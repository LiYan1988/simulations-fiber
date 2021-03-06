function param = derotate_constellation(param)
% Derotate the constellation back for visualization purpose

%% Declare some new variables
% Derotated points in the constellation
param.signal_received_constellation_derotate = ...
    cell(1, param.channel_number);
% Rotation angle from received to original constellations
param.angle_rotation = zeros(1, param.channel_number);
% Centers of point clouds after de-rotation
param.cloud_centers_derotation = cell(1, param.channel_number);

% % These variables will be used for now
% % Original constellation (after sorting) with [amplitude, angle]
% param.constellation_original = cell(1, param.channel_number);

% % Received constellation cloud centers (after sorting) with [amplitude,
% % angle]
% param.constellation_received = cell(1, param.channel_number);

for cidx=1:param.channel_number
    %% Calculate amplitude and angle of the original constellation
    
    % the symbols at transmitter side
    u = param.data_mod_symbol_channel{cidx};
    % Unique constellation points
    C_template = unique(u);
    % Amplitude of constellation points
    C_template_amplitude = sqrt(sum(abs(C_template).^2, 2));
    % Normalize amplitude of constellation points so that the maxium amplitude
    % equals 1
    C_template_amplitude = C_template_amplitude/max(C_template_amplitude);
    % Calculate angles of constellation points in [pi], i.e., the angles are
    % multiples of pi
    C_template_angle = angle(C_template)/pi;
    
    % Contains amplitudes and angles
    C_template = sortrows([C_template_amplitude, C_template_angle], [1, 2]);
    
    %% Real constellation cloud centers
    C_receive = param.cloud_centers{cidx};
    C_receive_amplitude = sqrt(sum(C_receive.^2, 2));
    C_receive_amplitude = C_receive_amplitude./max(C_receive_amplitude);
    
    % Remove noise in amplitude, this part can be automated by smarter
    % algorithm. Hard code for speed now.
    % This is needed because the received constellation cloud centers are
    % sorted by their amplitudes and angles
    if param.constellation_size(cidx) == 2
        a = C_receive_amplitude>0.5;
        b = C_receive_amplitude<0.5;
        C_receive_amplitude(a) = 1;
        C_receive_amplitude(b) = 0;
    elseif param.constellation_size(cidx) == 16
        a = C_receive_amplitude>0.8723;
        b = (C_receive_amplitude<0.8723)&(C_receive_amplitude>0.5395);
        c = C_receive_amplitude<0.5395;
        C_receive_amplitude(a) = 1;
        C_receive_amplitude(b) = 0.74;
        C_receive_amplitude(c) = 0.33;
    end
    
    C_receive_angle = angle(C_receive(:, 1)+1i*C_receive(:, 2))/pi;
    
    [C_receive, ~] = sortrows([C_receive_amplitude, C_receive_angle], [1, 2]);
    
    %% Calculate angle to de-rotation
    % to-do:
    % for ook, should calculate the relative rotation between the two point
    % clouds
    % for 16qam, should calculate the relative rotation between the point
    % clouds and the center of the constellation
    % but for now, it is okay to roughly rotate the constellation for
    % visualization purpose
    if param.constellation_size(cidx) == 2
        rotation_angle = mean(C_template(2, 2)-C_receive(2, 2));
    elseif param.constellation_size(cidx) == 16
        a = mean(C_template(1:4, 2)-C_receive(1:4, 2));
        b = mean(C_template(5:12, 2)-C_receive(5:12, 2));
        c = mean(C_template(13:end, 2)-C_receive(13:end, 2));
        rotation_angle = c;
    end
    
    param.angle_rotation(cidx) = rotation_angle;
    
    %% De-rotate received signal back in the constellation
    % the received signal, corresponds to points in the received constellation
    signal = param.signal_received_constellation{cidx};
    
    % convert to complex number
    signal_complex = signal(:, 1)+signal(:, 2)*1i;
    % de-rotate signal
    signal_complex = signal_complex*exp(1i*pi*rotation_angle);
    % convert back to Nx2 vector
    signal_derotate = zeros(size(signal));
    signal_derotate(:, 1) = real(signal_complex);
    signal_derotate(:, 2) = imag(signal_complex);
    
    param.signal_received_constellation_derotate{cidx} = signal_derotate;
    
    %% Plot de-rotated constellation
%     plot(signal_derotate(:, 1), signal_derotate(:, 2), '.')

    %% Find centers of cloud points
    % Could also apply de-rotation angle to cloud centers, which have been
    % calculated already
%     opts = statset('UseParallel', true);
%     [~, C] = kmeans(signal_derotate, param.constellation_size(cidx), ...
%         'Display', 'off', 'maxiter', 1000, ...
%         'Replicates', 64, 'Options', opts);
% 
%     param.cloud_centers_derotation{cidx} = C;
    
    u = param.cloud_centers{cidx};
    u = (u(:, 1)+1i*u(:, 2))*exp(1i*pi*rotation_angle);
    param.cloud_centers_derotation{cidx} = [real(u), imag(u)];
end