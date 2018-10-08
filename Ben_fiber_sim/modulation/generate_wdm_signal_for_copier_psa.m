function [p, u0]=generate_wdm_signal_for_copier_psa(p)
% Generates a WDM signal, and reserves space for idler waveforms (this
% assumes that a copier-PSA scheme is going to be used)

% power_chan is the pump power [dBm]

% Setting default values
if ~exist('p', 'var')
    p.exists = [];
    p = get_parameters(p);
end

write_log(p, now, 'Generating Signal (generate_wdm_signal_for_copier_psa)');
write_log(p, now, sprintf('Input CW power, %.2f dBm', p.light.power_dBm));

u0 = zeros(1, length(p.t)); % Pre-allocating output vector

% Generate and modulate light
u_cw = light_signal(p); % If not including laser phase noise, leave this outside of loop. Else, place inside loop
measure_power(u_cw, 'laser output');
p.tx.ideal = zeros(p.link.N_chan, p.N_symb);
for k=1:p.link.N_chan

    % Storing ideal transmitted constellation
    p.tx.ideal(k, :) = iq_modulator(p, ones(1, p.N_symb), [p.electric.levels(data_to_symnum(p.data(1, :, k)));p.electric.levels(data_to_symnum(p.data(2, :, k)))]);
    p.tx.ideal(k, :) = p.tx.ideal(k, :)./sqrt(mean(abs(p.tx.ideal(k, :)).^2)); % Normalizing
    
    % Set up the electric driving signals
    p.electric.data = p.data(1, :, 1, k); e_i = electric_signal(p);
    p.electric.data = p.data(2, :, 1, k); e_q = electric_signal(p);

    u_mod = iq_modulator(p, u_cw, [e_i; e_q]); % Modulating data onto a new channel TODO: Make sure that all of the data is independent!!!
    measure_power(u_mod, sprintf('after mod. (CH%d)', k));
    power_out = measure_power(u_mod);
        
    u_mod = u_mod.*exp(1j*p.w_index.sig(k)*p.t); % Shifting channel to new carrier frequency
    u_mod = set_power(u_mod, p.fopa.sig_in_dBm); % Setting signal power into FOPA
    
    u0 = u0 + u_mod; % Adding new channel to the output
    
    write_log(p, now, sprintf('Generated channel %d with %d symbols, Signal power after modulation, %.2f dBm, Center frequency, %.2f GHz (Relative to %.2f THz)', k, p.N_symb, power_out, p.w_index.sig(k)/(2*pi*1e9), p.const.c/p.const.lambda*1e-12));
end

% Convert to single precision; enables GPU acceleration, reduces disk-space
% of output files, but introduces an error that grows with the number of
% iterations in the split-step solver.
% Evaluate whether the extra precision is needed on a project-by-project
% basis.
% u0 = single(u0);

p = rmfield(p, 'electric'); % Cleaning up useless fields

write_log(p, now, 'Signal Generation Complete');