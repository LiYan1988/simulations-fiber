function plot_current_signal(param, scale)
% Plot signal in both time and frequency domains

move_mean_point = 1050;
move_mean_f = movmean(abs(param.data_mod_f_current).^2*1e12, move_mean_point);
move_mean_s = sprintf('movmean(PSD, %.1f GHz)', param.df/2/pi*move_mean_point/1e9);

figure;
subplot(2, 1, 1)
plot(param.t_plot, 1e3*abs(param.data_mod_t_current).^2)
xlabel('Time ($\mu$s)', 'interpreter', 'latex')
ylabel('Power (mW)', 'interpreter', 'latex')
grid on;
if strcmp(scale, 'log')
    subplot(2, 1, 2)
    hold on;
    h1 = plot(param.f_plot, 10*log10(abs(param.data_mod_f_current).^2*1e12), 'displayname', 'PSD');
    h2 = plot(param.f_plot, 10*log10(move_mean_f), 'displayname', move_mean_s, ...
        'linewidth', 2);
    xlabel('Frequency (GHz)', 'interpreter', 'latex')
    ylabel('PSD (dBm/GHz)', 'interpreter', 'latex')
    grid on;
    legend([h1, h2])
elseif strcmp(scale, 'linear')
    subplot(2, 1, 2)
    hold on;
    h1 = plot(param.f_plot, abs(param.data_mod_f_current).^2*1e12, 'displayname', 'PSD');
    h2 = plot(param.f_plot, move_mean_f, 'displayname', move_mean_s, ...
        'linewidth', 2);
    xlabel('Frequency (GHz)', 'interpreter', 'latex')
    ylabel('PSD (mW/GHz)', 'interpreter', 'latex')
    grid on;
    legend([h1, h2])
end
