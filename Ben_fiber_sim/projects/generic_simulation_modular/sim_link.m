function p = sim_link(slurm_id, varargin)
% Propagates a signal through a generic transmission link. Takes as inputs
% an identifier (a timestamp) for the transmitted signal and keyval pairs
% to set parameter sweeps. Saves the optical waveforms after every span in
% a results folder.
%
% INPUTS:
%    p : The parameter struct from caller
%
% Benjamin Foo, 2018-02-19
%
% Based on code by:
% Pontus Johannisson, 2010-04-19
% This software is distributed under the terms of the GNU General
% Public License version 2

%% Setting up paths
addpath('../../filter');
addpath('../../modulation');
addpath('../../nlse');
addpath('../../signal_plot');
addpath('../../transmission');
addpath('../../utils');

%% Loading location of parameter and transmitter files
try
    fid = fopen(sprintf('Job_%d_files.txt', slurm_id), 'r');
    sim_id = fscanf(fid, '%s');
    fclose(fid);
catch
    error('Invalid Job ID.')
end

%% Loading transmitter files
try
    % Loading parameter structure
    p   = load(sprintf('%s_Parameters.mat', sim_id), 'p');
    p   = p.p;
    
    if p.flag.b2b || (p.link.N_span == 0)
        return
    end
    
    p.t = 0:p.dt_samp:(p.N_samp - 1)*p.dt_samp;
    p.f = 1 /p.T_period*[(0:p.N_samp/2 - 1) (-p.N_samp/2:-1)];
        
    % Re-opening log
        if isfield(p.log, 'filename')
        p.log.fid = fopen(p.log.filename, 'a'); % Open Log file in append mode
            if p.log.fid == -1
                warning('Could not open log file')
            end
        end
catch
    if ~exist('sim_id', 'var')
        error('Please specify the ID of the transmitter files.')
    elseif ~exist(sprintf('../results/%s_Parameters.mat', sim_id), 'file')
        error('Could not find specified transmitter files. Please check that you have the correct ID.')
    else
        error('Unknown error...')
    end
end

% Dynamic parameter allocation using keyval pairs in the function call
% Particularly useful for parameter sweeps
p_check = struct('power', 0, 'convergence', 0, 'spans', 0, 'resume', 0, 'disp_map', 0); % Flags checking if a parameter is modified several times in the same function call. This should be adapted to include every case in the switch
p_str = '';
if nargin > 1
    for sweep_ind = 1:numel(varargin)
        try key = char(varargin(sweep_ind));
            switch lower(key)
                case 'power'
                    % Change WDM Launch power [dBm]
                    if ~p_check.power
                        value = cell2mat(varargin(sweep_ind+1));
                        if isnumeric(value)
                            p.link.wdm_pow = value;
                            fprintf('Setup: Launch power %d dBm\n', p.link.wdm_pow)
                            write_log(p, now, sprintf('Setup: Launch power %d dBm', p.link.wdm_pow)); % For convergence testing
                        end 
                        p_str = [p_str, sprintf('%d_dBm_', p.link.wdm_pow)];
                        p_check.power = 1; % Setting flag that power has been changed
                    else
                        warning('%s has previously been set in this function call. Additional attempts to set this parameter are ignored.', key)
                    end
                case 'convergence'
                    % Change number of steps per nonlinear length
                    if ~p_check.convergence
                        value = cell2mat(varargin(sweep_ind+1));
                        if isnumeric(value)
                            if value > 0
                                p.steps_per_L_NL = value;
                            else
                                p.steps_per_L_NL = 1;
                                warning('The number of steps per nonlinear length must be greater than 0. Setting to 1 instead.')
                            end
                            fprintf('Setup: Steps per nonlinear length %d\n', p.steps_per_L_NL)
                            write_log(p, now, sprintf('Setup: Steps per nonlinear length %d ', p.steps_per_L_NL)); % For convergence testing
                            p_str = sprintf('%s%d_nlse_steps_', p_str, p.steps_per_L_NL);
                        end
                        p_check.convergence = 1; % Setting flag that power has been changed
                    else
                        warning('%s has previously been set in this function call. Additional attempts to set this parameter are ignored.', key)
                    end
                case 'spans'
                    % Change number of steps per nonlinear length
                    if ~p_check.spans
                        value = cell2mat(varargin(sweep_ind+1));
                        if isnumeric(value)
                            if value > 0
                                p.link.N_span = value; 
                            else
                                p.link.N_span = 1; 
                                warning('The number of spans must be greater than 0. Setting to 1 instead.')
                            end
                            fprintf('Setup: Spans %d\n', p.link.N_span)
                            p.rx.edc_L = p.smf.L*p.link.N_span; % Remember to update any dependencies!
                            write_log(p, now, sprintf('Setup: Spans %d ', p.link.N_span)); % For convergence testing
                        end
                        p_check.spans = 1;
                    else
                        warning('%s has previously been set in this function call. Additional attempts to set this parameter are ignored.', key)
                    end
                case 'resume'
                    % Change number of steps per nonlinear length
                    if ~p_check.resume
                        value = cell2mat(varargin(sweep_ind+1));
                        if isnumeric(value)
                            if value > 0
                                start_span = value; 
                            else
                                start_span = 0; 
                                warning('Must start from a non-negative span index.')
                            end
                            fprintf('Setup: Resuming simulation from span %d\n', start_span)
                            write_log(p, now, sprintf('Setup: Resuming simulation from span %d ', start_span)); % For convergence testing
                        end
                        p_check.resume = 1;
                    else
                        warning('%s has previously been set in this function call. Additional attempts to set this parameter are ignored.', key)
                    end
                case 'dispersion'
                    % Change dispersion pre-compensation. Value is the
                    % equivalent amount of SMF in METERS(!) that should be
                    % pre-compensated.
                    if ~p_check.disp_map
                        value = cell2mat(varargin(sweep_ind+1));
                        if isnumeric(value)
                            p.psa.pre_comp_L = value;
                        end
                        p_str = sprintf('%s%.1f_km_precomp_', p_str, p.psa.pre_comp_L/1e3);
                        fprintf('Setup: Dispersion pre-compensation %.1f km\n',  p.psa.pre_comp_L/1e3)
                        write_log(p, now, sprintf('Setup: Dispersion pre-compensation %.2f km ', p.psa.pre_comp_L/1e3)); % For convergence testing
                        p_check.disp_map = 1;
                    else
                        warning('%s has previously been set in this function call. Additional attempts to set this parameter are ignored.', key)
                    end
                otherwise
                    warning('Unknown command ''%s''.', key);
            end
        catch
            % This is the value part. Do nothing.
        end
    end
end

% Some error checking
if ~exist('start_span', 'var')
    start_span = 0;
end
if start_span >= p.link.N_span
    % Trying to start a simulation which is not needed
    error('Total required link length occurs before the start point (spans <= resume)')
end

% Loading input waveform
if ~start_span % start_span = 0 implies load from TX
    load(sprintf('%s_Modulator_output_optical_waveform.mat', sim_id), 'u0')
    
    % Adding some noise to mimic transceiver impairments
    %     if isfield(p.tx, 'snr_db')
    %         u0 = awgn(u0, p.tx.snr_db, 'measured');
    %     end

    % Set signal input power assuming that all the channels have the same input power.
    if isfield(p, 'wdm_power')
        u = set_power(u0, p.link.wdm_pow + 10*log10(p.link.N_chan));
        measure_power(u, 'WDM input power');
    else
    u=u0;
    end
    
    % Adding laser phase noise
    if isfield(p.tx, 'linewidth')
        u = add_phase_drift(u, 2*pi*p.tx.linewidth*p.dt_symb);
    end
    
    if p.flag.osa % Monitoring transmitted spectrum after power set
        fig_handle_spec = figure(81); clf(fig_handle_spec); plot_spectrum(p, u); title(fig_handle_spec.CurrentAxes, 'Fiber Input Optical Spectrum');
        fig_handle_spec.Position = [2560 200 700 500];
    end
    
else
    try
        load(sprintf('%s_Fiber_optical_waveform_%sspan_%d.mat', sim_id, p_str, start_span), 'u')
    catch
        error('Could not find file: %s.', sprintf('%s_Fiber_optical_waveform_%sspan_%d.mat', sim_id, p_str, start_span))
    end
end

fprintf('%s, Link started: %s %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), datestr(p.timestamp, 'yyyy-mm-dd_HH.MM.SS'), mfilename);
write_log(p, p.timestamp, sprintf('%s Started', mfilename));

%% Propagate through a fiber
measure_power(u, 'into link.  ');
p.link.wdm_pow=measure_power(u);
if isfield(p, 'wdm_power')
    write_log(p, now, sprintf('Target per-channel launch power is %d dBm', p.link.wdm_pow));
else
    write_log(p, now, sprintf('WDM power is %d dBm', p.link.wdm_pow));
end

% if p.flag.psa
%     % Copier emulation
%     u_idler = conj(u);
%     
%     % Dispersion pre-compensation
%     u = dispersive_propagation(u, p.omega, -p.psa.pre_comp_L*p.smf.D);
%     u_idler = dispersive_propagation(u_idler, p.omega, -p.psa.pre_comp_L*p.smf.D);
% end

%Calculating the number of delta Z to solve the NLSE, put in the loop if you want it to change after every span.
%p.smf.delta_z = 1/max(max(abs(u).^2))/p.smf.gamma/p.steps_per_L_NL;
p.smf.delta_z = p.smf.L/p.steps_per_L_NL; % Different implementation to set solver step-size

% Transmit through fiber spans
for span_ind = start_span+1:p.link.N_span
    rng(mod(p.rng.seed + span_ind, 2^32-1)); % Re-setting RNG every span
   
    p.fiber = p.smf;
    if isfield(p.flag, 'gpu') && p.flag.gpu
        try
            u = nlse_spm_gpu(p, u);
        catch
            warning('Could not find compatible GPU. Running solver without GPU acceleration')
            write_log(p, now, 'You forgot to disable the GPU flag when running this on the cluster (again)... idiot. Running CPU-only solver')
            u = nlse_spm(p, u); % Running CPU-only NLSE solver
            p.flag.gpu = 0; % Setting flag to 0 so that this doesn't happen again
        end
    else
        %fprintf('Iterations, %d steps per L NL,%d\n',p.smf.L/p.steps_per_L_NL,p.steps_per_L_NL); %If it changes all the time is good to seehow it changes
        u = nlse_spm(p, u);
    end
    
    % Propagating idler
%     if p.flag.psa
%         if isfield(p.flag, 'gpu') && p.flag.gpu
%             u_idler = nlse_spm_gpu(p, u_idler);
%         else
%             u_idler = nlse_spm(p, u_idler);
%         end
%     end
    
    measure_power(u, sprintf('after span %d', span_ind));

%     if p.flag.psa
%         % Dispersion post-compensation
%         u = dispersive_propagation(u, p.omega, -(p.smf.L-p.psa.pre_comp_L)*p.smf.D);
%         u_idler = dispersive_propagation(u_idler, p.omega, -(p.smf.L-p.psa.pre_comp_L)*p.smf.D);
%         p.rx.edc_L = 0;
%     end
    
    % Saving fiber output
    if ~mod(span_ind, p.file.save_span_incr)
        if p.flag.psa
            save(sprintf('%s%s_Fiber_optical_waveform_%sspan_%d.mat', p.file.path, datestr(p.timestamp, 'yyyy-mm-dd_HH.MM.SS'), p_str, span_ind), 'u', 'u_idler')
        else
            save(sprintf('%s%s_Fiber_optical_waveform_%sspan_%d.mat', p.file.path, datestr(p.timestamp, 'yyyy-mm-dd_HH.MM.SS'), p_str, span_ind), 'u')
        end
    end
       
    write_log(p, now, sprintf('SMF: Propagated %d km', span_ind*p.smf.L/1e3));
    
    if span_ind < p.link.N_span
        u = edfa(u, p.fiber.L*p.fiber.att_db, p.link.edfa_nf_db, p.const.nu, p.dt_samp); %Amplifying after every span.
        if p.flag.edc_es==1
            u = dispersive_propagation(u,p.omega,-p.smf.D*p.smf.L); %Dispersion compensation at each span end.
        end
    end
end % for span_ind...

write_log(p, now, sprintf('Transmission through %d km SMF completed', p.link.N_span*p.smf.L/1e3));
% 
% fig_handle_spec_out = figure(82); clf; plot_spectrum(p, u); title(fig_handle_spec_out.CurrentAxes, 'Output Spectrum'); % Output Spectrum
% fig_handle_spec_out.Position = [3260 850 700 500];

p = rmfield(p, {'fiber', 'f', 't'}); % Removing useless fields

fprintf('%s Link finished: %s %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), datestr(p.timestamp, 'yyyy-mm-dd_HH.MM.SS'), mfilename);
write_log(p, now, 'Dispersion Management Simulation Complete\n');
if p.log.fid ~= -1 % If log is open, close it
    fclose(p.log.fid);
end

%% Removing paths
rmpath('../../filter');
rmpath('../../modulation');
rmpath('../../nlse');
rmpath('../../signal_plot');
rmpath('../../transmission');
rmpath('../../utils');
