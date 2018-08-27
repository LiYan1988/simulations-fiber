clc;
clear;
close all;

% Write functions to implement:
% 1. find centers of point clouds in the constellation diagram
% 2. de-rotate the constellation disgram back for visualization purpose

%% Load data
load ssf_signal_polarization_5.mat

%% Test function
% param = center_constellation(param);

%% Signal clusters
close all;
clc;

cidx = (N-1)/2+1;
[xt_dc, ~, ~] = dispersion_compensation(param.data_mod_t_current, cidx,...
    param.beta2, param.beta3, 0, param);
xt_dc = xt_dc(param.delay_filter_channel(cidx)+1:...
    end-param.delay_filter_channel(cidx));

signal = zeros(size(xt_dc, 1), 2);
signal(:, 1) = real(xt_dc);
signal(:, 2) = imag(xt_dc);

signal = downsample(signal, param.sample_per_symbol(cidx), ...
    param.shift_channel_time(cidx));

% Find center of points 
opts = statset('UseParallel', true);
[idx,C,sumd,D] = kmeans(signal, 16, 'Display', 'final', ...
    'maxiter', 1000, 'Replicates', 64, 'Options', opts);

figure;
hold on;
box on;
grid on;
plot(signal(:, 1), signal(:, 2), '.')
plot(C(:, 1), C(:, 2), 'x')

centers = C(idx, :);
noise = sqrt(sum((signal - centers).^2, 2));

noise_power = zeros(16, 1);
for k=1:16
    noise_power(k) = mean(noise(idx==k));
end

snr_point_cloud = sqrt(sum(C.^2, 2))./noise_power;
snr_total = mean(sqrt(sum(centers.^2, 2)))/mean(noise);

%% 
cidx = 6;
u = param.data_mod_symbol_channel{cidx};
% scatterplot(u)

C_template = unique(u);
C_template_amplitude = sqrt(sum(abs(C_template).^2, 2));
C_template_amplitude = C_template_amplitude/max(C_template_amplitude);
C_template_angle = angle(C_template)/pi;

C_template = sortrows([C_template_amplitude, C_template_angle], [1, 2]);

%% Re-rotate the constellation diagram
C_amplitude = sqrt(sum(C.^2, 2));
C_amplitude = C_amplitude./max(C_amplitude);

C_amplitude(C_amplitude>0.8723) = 1;
C_amplitude((C_amplitude<0.8723)&(C_amplitude>0.5395)) = 0.74;
C_amplitude(C_amplitude<0.5395) = 0.33;

C_angle = angle(C(:, 1)+1i*C(:, 2))/pi;

[C_constellation, C_idx] = sortrows([C_amplitude, C_angle], [1, 2]);
C = C(C_idx, :);

% figure;
% hold on;
% plot(C_constellation(:, 1))
% plot(C_constellation(:, 2))
% 
% figure;
% plot(C_constellation(:, 2)-C_template(:, 2))
% title('Different constellation points have different ratation angles')

rotation_inner = mean(C_constellation(1:4, 2));
rotation_midium = mean(C_constellation(5:12, 2));
rotation_outer = mean(C_constellation(13:end, 2));

fprintf('Inner rotation is %.2fpi\n', rotation_inner)
fprintf('Medium rotation is %.2fpi\n', rotation_midium)
fprintf('Outer rotation is %.2fpi\n', rotation_outer)

%% Retate constellation back
% use rotation_medium to rotate the constellation back
signal_complex = signal(:, 1)+signal(:, 2)*1i;
signal_complex = signal_complex*exp(-1i*pi*rotation_midium);

signal_derotate = zeros(size(signal));
signal_derotate(:, 1) = real(signal_complex);
signal_derotate(:, 2) = imag(signal_complex);

% plot the constellation again
figure;
hold on;
plot(signal_derotate(:, 1), signal_derotate(:, 2), '.')
% plot(real(u), imag(u), 'x')