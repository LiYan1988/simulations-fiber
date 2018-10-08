function [p] = coherent_receiver(p, u)
% Implementation of a coherent receiver to recover a signal
%
% 'Electronic' CD compensation
% Re-samples to [p.rx.os] times oversampling 
% Dynamic Equalizer
% EVM calculation
% BER calculation
% GMI calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Input monitoring stage
write_log(p, now, 'Started coherent receiver');

% Resetting RNG seed so that the same thermal noise is generated for all
% distances (makes a fairer comparison) - should disable this if a
% statistical test is preferred (i.e. different noise each time the script is called)
rng(mod(p.rng.seed - 1, 2^32-1));

% Monitoring received optical spectrum
rx_power = measure_power(u);
measure_power(u, 'received signal power ');
write_log(p, now, sprintf('RX: Received signal power %.2f dBm', rx_power));

if p.flag.osa % Monitor receive optical spectrum
    rx_power = ceil(rx_power/10)*10;
    fig_handle = figure(91); clf(fig_handle); plot_spectrum(p, u); axis(fig_handle.CurrentAxes, [-p.rx.filter_bandwidth_hz/1e9, p.rx.filter_bandwidth_hz/1e9, rx_power-70, rx_power]); title(fig_handle.CurrentAxes, 'Received Optical Spectrum')
    fig_handle.Position = [2060, 40, 560, 420];
end


%% Add thermal noise in receiver
if isfield(p.rx, 'thermal_noise') && isfield(p.rx, 'responsivity')
    %     % Add AWGN, the complex noise amplitude should be such that the noise
%     % power is P_n
%     u_n = sqrt(P_n)/sqrt(2)*(randn(size(u)) + 1i*randn(size(u))); % From
%     EDFA function

    noise = p.rx.thermal_noise*sqrt(p.f_samp/2)*(randn(1, length(u)) + 1i*randn(1, length(u)));
    demux = struct('delta_f', p.rx.filter_bandwidth_hz, 'f0', 0, 'f', p.f);
    noise = filter_ideal_bp(demux, noise); % Only adds noise in the receiver's electrical bandwidth. TODO: Consider using a different type of filter to emulate the photodiode?
    u = p.rx.responsivity*u + noise; 

    % Monitoring received spectrum after thermal noise loading
    measure_power(filter_ideal_bp(demux, noise), 'thermal noise ');
    measure_power(filter_ideal_bp(demux, u), 'total received power ');
end
if p.flag.osa % Monitor received electrical spectrum (after adding thermal noise)
    fig_handle = figure(92); clf(fig_handle); plot_spectrum(p, u); pause(0.5); axis(fig_handle.CurrentAxes, [-p.rx.filter_bandwidth_hz/1e9, p.rx.filter_bandwidth_hz/1e9, rx_power-70, rx_power]); title(fig_handle.CurrentAxes, 'Received Electrical Spectrum (with thermal noise)')
    fig_handle.Position = [2620, 40, 560, 420]; 
end

% noise_power = measure_power(filter_ideal_bp(p, noise));
% write_log(p, now, sprintf('RX: Thermal noise %.2f dBm', noise_power));

%% Emulating electronic dispersion compensation with dispersive fiber (optional)
if (isfield(p.rx, 'edc_L')) && (p.rx.edc_L ~= 0) && isfield(p.flag, 'b2b') && ~p.flag.b2b
    % Run EDC if a) the parameter edc_L has been specified as non-zero (can
    % be negative), and b) this is not a back-to-back simulation.
    
    write_log(p, now, sprintf('RX: Applying EDC for %.2f km', p.rx.edc_L/1e3));
    u = dispersive_propagation(u,p.omega,-p.smf.D*p.rx.edc_L);
    
    % Mean group delay compensation
    u = circshift(u, round(-p.link.N_span*p.smf.L*p.smf.D*p.w_index.sig(p.rx.chan)/(2*pi)*p.const.lambda^2/p.const.c/p.dt_samp));
end

%% Downsampling with ADCs and dynamic equalizer
% DD-LMS equalizer (from Attila's Experimental code)
if isfield(p, 'rx') && isfield(p.rx, 'equalizer') % Try and run a dynamic equalizer
    try
        % Sampling noisy electrical signal (ADCs)
        p.rx.u = resample(double(u), p.rx.os, p.samp_per_symb); % Sampling to desired rate for DSP NOTE: Apparently, resample only works on doubles
        p.rx.u = p.rx.u/sqrt(mean(abs(p.rx.u).^2));
        
        if isfield(p.rx, 'dc_bias') && isfield(p.rx.equalizer,'Iterations_DDLMS') && (p.rx.equalizer.Iterations_DDLMS > 0)
            p.rx.u = p.rx.u + p.rx.dc_bias; % Adding a small carrier for Frequency offset estimation
        end
                
        write_log(p, now, 'RX: DDLMS Equalizer');
        write_log(p, now, sprintf('RX: Samples per symbol %d', p.rx.os));
        write_log(p, now, sprintf('RX: Equalizer taps %d', p.rx.equalizer.NTap));
        write_log(p, now, sprintf('RX: CMA iterations %d \t RD-CMA iterations %d \t DD-LMS iterations %d', p.rx.equalizer.Iterations_CMA, p.rx.equalizer.Iterations_RDCMA, p.rx.equalizer.Iterations_DDLMS));
        write_log(p, now, sprintf('RX: BPS test angles %d', p.rx.os));

        [p.rx.u, ~, p.rx.taps, ~, ~, ~, p.rx.f_offset, ~, p.rx.phase] = DDLMS_eq(p.rx.equalizer, single(p.rx.u.')); % DD-LMS algorithm expects input to be a single-precision column vector NOTE: DD-LMS algorithm removes some symbols from front and end of data stream. May need to re-synchronize outputs
        p.rx.u = p.rx.u.';

        write_log(p, now, sprintf('RX: Estimated Frequency Offset %.4f', p.rx.f_offset));
    catch
        warning('Incorrect inputs to DD-LMS equalizer code. Reverting to simple sampling case')
        error_flag = 1;
        
        % If not running an equalizer, downsample to one sample per symbol
        p.rx.u = u(p.idx_symb); % 'Ideal' downsampling to one sample per symbol
        p.rx.u = p.rx.u/sqrt(mean(abs(p.rx.u).^2));
    end
else
    
    p.rx.u = u(p.idx_symb);
    p.rx.u = p.rx.u/sqrt(mean(abs(p.rx.u).^2));
end
% 
% % Use external phase recovey (4th-power) code if 1) no equalizer is specified, 2) DD-LMS has 0 iterations or 3) the equalizer was not run because of an error
if ~isfield(p.rx, 'equalizer') || (isfield(p.rx, 'equalizer') && isfield(p.rx.equalizer, 'Iterations_DDLMS') && ~p.rx.equalizer.Iterations_DDLMS) || (exist('error_flag', 'var') && error_flag)   
    write_log(p, now, 'RX: Viterbi-Viterbi Phase Recovery');
    [p.rx.u, ~] = CPE_4thPower_SlidingWindow([p.rx.u(end-p.rx.vv_block_length/2+1:end), p.rx.u, p.rx.u(1:p.rx.vv_block_length/2)], 0, p.rx.vv_block_length); % Adding samples at front and end of received symbols (assuming periodicity) to negate symbol removal from CPE TODO: Kind of works, but may want a better solution eventually.
    p.rx.u = p.rx.u.';
end

%% Calculating Received Signal Quality (BER, EVM, GMI)
write_log(p, now, 'RX: Calculating Signal Performance Metrics');
% Testing pi/2 phase rotations to minimize BER
test_ber = ones(1, 4);
for test_ind = 1:4
    test_u = p.rx.u.*exp(1j*pi/2*test_ind);
    test_data = bit_demapper(p, test_u);
    
    test_data = reshape(test_data, 1, p.bits_per_symbol*length(p.rx.u));
    tx_data = reshape(p.data(:, :, 1, p.rx.chan), 1, p.bits_per_symbol*p.N_symb); % TODO: Make polarization compatible; Make WDM compatible
    
    test_ber(test_ind) = ber_counter(test_data, tx_data);
end

[p.rx.ber, ber_ind] = min(test_ber);
p.rx.u = p.rx.u.*exp(1j*pi/2*ber_ind);

% Storing binary data
p.rx.data = bit_demapper(p, p.rx.u);

% % Syncing constellations for EVM and GMI calculations
[val, lags] = crosscorr(real(p.rx.u), real(p.tx.ideal(p.rx.chan, :)), length(p.rx.u)-1);
[~, sync_ind] = max(abs(val));
sync = lags(sync_ind);
tx_sync = circshift(p.tx.ideal(p.rx.chan, :), -sync);

% Calculating EVM
error = p.rx.u - tx_sync(1:length(p.rx.u));
p.rx.evm_db = 10*log10(sqrt(mean(abs(error).^2))/sqrt(mean(abs( tx_sync(1:length(p.rx.u)) ).^2)));

% GMI calculation
p.rx.gmi = calcGMI_withNormalization(tx_sync(1:length(p.rx.u)), p.rx.u);

fprintf('RX: BER = %.2e \t EVM = %.2f dB \t GMI = %.2e bit/Symb\n', p.rx.ber, p.rx.evm_db, p.rx.gmi);

%% Saving results
result = struct('ber', p.rx.ber, 'evm_db', p.rx.evm_db, 'gmi', p.rx.gmi, 'u0', p.tx.ideal, 'u', p.rx.u);
% save_file_name = [datestr(p.timestamp, 'yyyy-mm-dd_HH.MM.SS'), '_Simulation_result_channel_', num2str(p.rx.chan, '%d'), '_sweep_iteration', num2str(round(p.sweep.ind)), '.mat'];
% save(['../results/', save_file_name], 'result');

write_log(p, now, sprintf('RX: BER %.2e', result.ber));
write_log(p, now, sprintf('RX: EVM %.2f dB', result.evm_db));
write_log(p, now, sprintf('RX: GMI %.4f bits/symbol', result.gmi));

if p.flag.rx_perf
    fig_handle = figure(100); clf(fig_handle); fig_handle.CurrentAxes = axes; plot(fig_handle.CurrentAxes, real(p.rx.u), imag(p.rx.u), '.'); axis(fig_handle.CurrentAxes, [-1.5, 1.5, -1.5, 1.5]); grid(fig_handle.CurrentAxes, 'ON');
    title(fig_handle.CurrentAxes, 'RX Constellation');
    xlabel(fig_handle.CurrentAxes, ['BER = ', num2str(p.rx.ber, '%.2e'), '  EVM = ', num2str(p.rx.evm_db, '%.2f'), ' dB   GMI = ', num2str(p.rx.gmi, '%.2e'), ' bit/Symb']);
    fig_handle.Position = [3180, 40, 560, 420]; 
end

write_log(p, now, 'Coherent receiver complete');
end

function [output_bits] = bit_demapper(p, u)
% Demaps the symbols to bits based on the iq_modulator function
% This is super crude. I should fix it when I have time.
    i_rx = real(u);
    q_rx = imag(u);
    
    i_lvl = round((p.bits_per_symbol-1)/2*(i_rx-min(i_rx)));
    q_lvl = round((p.bits_per_symbol-1)/2*(q_rx-min(q_rx)));
    
    i_bits = zeros(ceil(p.bits_per_symbol/2), length(u));
    q_bits = zeros(floor(p.bits_per_symbol/2), length(u));
    switch p.modulation.format
        case '16-QAM'
            for ind = 1:length(u)
                if i_lvl(ind) == 3
                    i_bits(:, ind) = [1; 1];
                elseif i_lvl(ind) == 2
                    i_bits(:, ind) = [0; 1];
                elseif i_lvl(ind) == 1
                    i_bits(:, ind) = [1; 0];
                else % i_lvl == 0
                    i_bits(:, ind) = [0; 0];
                end
                
                if q_lvl(ind) == 3
                    q_bits(:, ind) = [1; 1];
                elseif q_lvl(ind) == 2
                    q_bits(:, ind) = [0; 1];
                elseif q_lvl(ind) == 1
                    q_bits(:, ind) = [1; 0];
                else % l_lvl == 0
                    q_bits(:, ind) = [0; 0];
                end
                            
            end
        case 'QPSK'
            % Demapping bits (QPSK)
            for ind = 1:length(u)
                if i_rx(ind) > 0
                   i_bits(ind) = 0;
                else
                    i_bits(ind) = 1;
                end

                if q_rx(ind) > 0
                    q_bits(ind) = 0;
                else
                    q_bits(ind) = 1;
                end
            end
        case 'OOK'
            i_bits = i_lvl;
        case 'BPSK'
            i_bits = i_lvl;
        otherwise
            warning('Not Found!')
    end
    output_bits = [i_bits;q_bits];
end