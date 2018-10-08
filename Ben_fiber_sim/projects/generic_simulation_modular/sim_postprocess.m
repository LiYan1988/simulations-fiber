function q = sim_postprocess(file_path, sim_id, rx_chan, psa_flag)
% Post-processes the results from a cluster simulation.
% DESIGNED FOR USE ON LOCAL DESKTOP!
% This function attempts to look for ALL output files belonging to a
% simulation, and then searches for keywords in the filename to determine
% parameters that are important for post-processing (e.g. number of spans,
% amount of in-line optical dispersion compensation). 
%
% In this case, the post-processing emulates a phase-sensitive
% pre-amplifier by performing coherent superposition followed by an EDFA
% with a noise figure that is lower than the regular EDFA by 3-dB.
%
% INPUTS:
%    p : The parameter struct from caller
%
% Benjamin Foo, 2018-02-21

%% Setting up paths
addpath('../../filter');
addpath('../../modulation');
addpath('../../nlse');
addpath('../../signal_plot');
addpath('../../transmission');
addpath('../../utils');

%% Loading transmitter files
try
    % Loading parameter structure
    p   = load(sprintf('%s%s_Parameters.mat', file_path, sim_id), 'p');
    p   = p.p;
    
    if p.flag.b2b || (p.link.N_span == 0)
        q = [];
        return
    end
    p.rx.chan = rx_chan;
    
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
    elseif ~exist(sprintf('%s%s_Parameters.mat', file_path, sim_id), 'file')
        error('Could not find specified transmitter files. Please check that you have the correct ID.')
    else
        error('Unknown error...')
    end
end

% Checking for PIA or PSA operation for processing. Simulation should (almost)ALWAYS be run in PSA mode, and then post-processing should 
% perform coherent superposition.
% if ~exist('psa_flag','var')
%     psa_flag = 0;
%     func_stack = dbstack;
%     warning('Called from ''%s'': User did not specify whether or not to process a PSA result.', func_stack(2).name)
% end

fprintf('%s, Post-processing started: %s %s\n', datestr(p.timestamp, 'yyyy-mm-dd HH:MM:SS'), datestr(p.timestamp, 'yyyy-mm-dd_HH.MM.SS'), mfilename);
write_log(p, p.timestamp, sprintf('%s Started', mfilename));

% Obtaining list of RX waveforms
file_list = dir(sprintf('%s%s_Fiber_optical_waveform_*.mat', file_path, sim_id));

for file_ind = 1:numel(file_list)
    if psa_flag
         load(sprintf('%s/%s', file_list(file_ind).folder, file_list(file_ind).name), 'u', 'u_idler')
    else
        load(sprintf('%s/%s', file_list(file_ind).folder, file_list(file_ind).name), 'u')
    end
    
    fprintf('Processing %s\n', file_list(file_ind).name)
    write_log(p, now, sprintf('Processing %s', file_list(file_ind).name));
    
    seperator_ind = strfind(file_list(file_ind).name, '_'); % Finds all instances of underscore (used to separate inputs);
    seperator_ind = seperator_ind(5:end); % Removes the first 4 underscores - 1 used in sim ID, '_Fiber_optical_waveform' gives 3 more
                                          %                                                      ^     ^       ^
    seperator_ind = [seperator_ind, length(file_list(file_ind).name)-3]; % Adding the index for the start of '.mat'
    
    % Finding the number of spans
    str_start = strfind(file_list(file_ind).name, 'span');
    tmp_ind = find(seperator_ind>str_start); % Finds indices for seperators after the target string
    str = file_list(file_ind).name(seperator_ind(tmp_ind(1))+1:seperator_ind(tmp_ind(2))-1); % Takes the string between the next two seperators
    p.link.N_span = str2double(str);
    
    if isfield(p.flag, 'psa') && p.flag.psa
        p.rx.edc_L = 0;
    else
        p.rx.edc_L = p.link.N_span*p.smf.L;
    end
    
    try
        % Finding per-channel power
        str_start = strfind(file_list(file_ind).name, 'dBm');
        tmp_ind = find(seperator_ind<str_start); % Finds indices for seperators before the target string
        str = file_list(file_ind).name(seperator_ind(tmp_ind(end-1))+1:seperator_ind(tmp_ind(end))-1); % Takes the string between the next two seperators
        p.link.wdm_pow = str2double(str);
    catch
    end
    
    try
        % Finding steps per nonlinear length
        str_start = strfind(file_list(file_ind).name, 'nlse');
        tmp_ind = find(seperator_ind<str_start); % Finds indices for seperators before the target string
        str = file_list(file_ind).name(seperator_ind(tmp_ind(end-1))+1:seperator_ind(tmp_ind(end))-1); % Takes the string between the next two seperators
        p.steps_per_L_NL = str2double(str);
    catch
    end
    
    try
        % Finding dispersion pre-compensation
        str_start = strfind(file_list(file_ind).name, 'km_precomp');
        tmp_ind = find(seperator_ind<str_start); % Finds indices for seperators before the target string
        str = file_list(file_ind).name(seperator_ind(tmp_ind(end-1))+1:seperator_ind(tmp_ind(end))-1); % Takes the string between the next two seperators
        p.psa.pre_comp_L = str2double(str)*1e3;
    catch
    end
    
    clear seperator_ind
    
    %% Pre-amplifier
    if psa_flag
        [u, u_idler] = coherent_superposition(u, u_idler);
        u = edfa(u, p.smf.L*p.smf.att_db, p.link.edfa_nf_db-3, p.const.c/p.const.lambda, p.dt_samp); % Assumes the PSA has a 3-dB noise figure improvement
    else
        u = edfa(u, p.smf.L*p.smf.att_db, p.link.edfa_nf_db, p.const.c/p.const.lambda, p.dt_samp); 
    end

    
    %% Coherent receiver
    fig_handle_spec_out = figure(82); clf; plot_spectrum(p, u); title(fig_handle_spec_out.CurrentAxes, 'Coherent Receiver Input Spectrum'); % Output Spectrum
    fig_handle_spec_out.Position = [3260 850 700 500];
    
    p_tmp = coherent_receiver(p, u);
    
    % Saving a result struct (struct of arrays, for direct use in plotting)
    q.power(file_ind) = p_tmp.link.wdm_pow;
    q.span(file_ind) = p_tmp.link.N_span;
    q.nlse_steps(file_ind) = p_tmp.steps_per_L_NL;
    q.precomp(file_ind) = p_tmp.psa.pre_comp_L/1e3;
    q.chan(file_ind) = p_tmp.rx.chan;
    q.ber(file_ind) = p_tmp.rx.ber;
    q.evm_db(file_ind) = p_tmp.rx.evm_db;
    q.gmi(file_ind) = p_tmp.rx.gmi;
    

    
    % Saving a second result struct (an array of structs), to be saved in a
    % file in case further analysis is required
    temp = struct('power', p_tmp.link.wdm_pow, 'span', p_tmp.link.N_span, 'chan', p_tmp.rx.chan, 'nlse_steps', p_tmp.steps_per_L_NL, 'ber', p_tmp.rx.ber, 'evm_db', p_tmp.rx.evm_db, 'gmi', p_tmp.rx.gmi);
    if p.flag.rx_perf
%         results(file_ind) = struct('power', p_tmp.link.wdm_pow, 'span', p_tmp.link.N_span, 'chan', p_tmp.rx.chan, 'nlse_steps', p_tmp.steps_per_L_NL, 'ber', p_tmp.rx.ber, 'evm_db', p_tmp.rx.evm_db, 'gmi', p_tmp.rx.gmi, 'constellation', p_tmp.rx.u, 'data', p_tmp.rx.data);
        temp.constellation = p_tmp.rx.u;
        temp.data = p_tmp.rx.data;
        q.constellation(file_ind) = {p_tmp.rx.u};
        q.data(file_ind) = {p_tmp.rx.data};
    end
    results(file_ind) = temp;
end

try
    span_elements = numel(unique(q.span));
    power_elements = numel(unique(q.power));
    fields = fieldnames(q);

    % Reshaping array; Each column is one launch power, each row is a given
    % distance
    for field_ind = 1:numel(fields)
        q.(char(fields(field_ind))) = reshape(q.(char(fields(field_ind))), span_elements, power_elements);
    end
catch
    warning('Could not reshape results array. Take care with processing results')
end
% Save parameter and results structs
p = rmfield(p, {'t', 'f', 'omega'}); % Removing space-hungry fields
save(sprintf('%s%s_Results.mat', file_path, sim_id), 'p', 'results');

fprintf('%s RX finished: %s %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), datestr(p.timestamp, 'yyyy-mm-dd_HH.MM.SS'), mfilename);
write_log(p, now, 'Post-processing Complete\n');
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
