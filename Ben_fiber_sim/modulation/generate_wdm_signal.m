function [p, u0]=generate_wdm_signal(p)
% Generates a WDM signal

% Setting default values
if ~exist('p', 'var')
    p.exists = [];
    p = get_parameters(p);
end

u0 = zeros(1, length(p.t)); % Pre-allocating output vector
p.tx.ideal = zeros(p.link.N_chan, p.N_symb);
u_cw = light_signal_0dBm(p); %Linear value. If not including laser phase noise, leave this outside of loop. Else, place inside loop. %0dBm is always setted in light_signal

% Generate electrical signal and modulate light

for k=1:p.link.N_chan
    write_log(p, now, 'Generating Signal (generate_wdm_signal)');
    write_log(p, now, sprintf('Input CW power, %.2f dBm', p.power_array(k)));
    if p.link.modulations(k)==1
        %Defining signal modulation 
        p.modulation.format = 'OOK'; p.bits_per_symbol = 1;
        p.data(k).bits = round(rand(p.bits_per_symbol, p.N_symb, (p.flag.dual_pol>0)+1));
        e_levels       = [0.5, 0];
        p.electric.levels = e_levels;
        %Electrical Filter
        e_filter.type      = 'gaussian'; %ideal_bp,gaussian.
        e_filter.bandwidth = 0.65*p.f_symb(k);  % Filter bandwidth [Hz]
        e_filter.bw_spec   = 'HWHM';         % Half-width, half maximum
        p.electric.filter = e_filter;
   
    elseif p.link.modulations(k)==3
        %Defining signal modulation 
        p.modulation.format = 'QPSK'; p.bits_per_symbol = 2;
        p.data(k).bits = round(rand(p.bits_per_symbol, p.N_symb, (p.flag.dual_pol>0)+1));
        e_levels       = [0, 1];
        p.electric.levels = e_levels;
        
        %Electrical Filter
        e_filter.type      = 'gaussian';  %ideal_bp,gaussian.
        e_filter.bandwidth = 0.65*p.f_symb(k);  % Filter bandwidth [Hz]
        e_filter.bw_spec   = 'HWHM';         % Half-width, half maximum
        p.electric.filter = e_filter;
        
        %Equalizer Parameters for Phase recovery
        p.rx.equalizer.Iterations_DDLMS = 0; % Number of DD-LMS iterations
        p.rx.equalizer.StepSize_DDLMS   = 1e-6; % Step size for DD-LMS equalizer % NOTE: Step size can either be the same for every DDLMS iteration (single number) or different (vector with a value for the step size for each iteration)
        p.Constellation = [-1-1j 1-1j 1+1j -1+1j];
        if numel(p.rx.equalizer.StepSize_DDLMS) ~= 1 && numel(p.rx.equalizer.StepSize_DDLMS) ~= p.rx.equalizer.Iterations_DDLMS
            error('The DDLMS step size is incorrect. Please specify one step size for all iterations or a vector for the step size of each iteration.')
        end
        [~,ix]=sort(angle(p.Constellation));
        p.Constellation=p.Constellation(ix);
        p.rx.equalizer.Constellation= single(p.Constellation/sqrt(mean(abs(p.Constellation).^2)));
        
    elseif p.link.modulations(k)==4
        %Defining signal modulation 
        p.modulation.format = '16-QAM'; p.bits_per_symbol = 4;
        p.data(k).bits = round(rand(p.bits_per_symbol, p.N_symb, (p.flag.dual_pol>0)+1));
        
        %I'm using the Communications toolbox in order to acheive the signal with RRC filtering.
        %e_levels       = [0 acos(1/3)/pi acos(-1/3)/pi 1]; % Modulation levels
        %p.electric.levels = e_levels;
%       e_filter.type      = 'gaussian';     %ideal_bp,gaussian.
%       e_filter.bandwidth = 0.65*p.f_symb(k);  % Filter bandwidth [Hz]
%       e_filter.bw_spec   = 'HWHM';         % Half-width, half maximum
%       p.electric.filter = e_filter;
        
        %Equalizer Parameters for Phase recovery
        p.rx.equalizer.Iterations_DDLMS = 5; % Number of DD-LMS iterations
        p.rx.equalizer.StepSize_DDLMS   = [1e-3, 1e-3, 1e-4, 1e-5, 1e-6]; % Step size for DD-LMS equalizer % NOTE: Step size can either be the same for every DDLMS iteration (single number) or different (vector with a value for the step size for each iteration)
        p.Constellation=[-3-3i,-3-1i,-3+1i,-3+3i,-1-3i,-1-1i,-1+1i,-1+3i,1-3i,1-1i,1+1i,1+3i,3-3i,3-1i,3+1i,3+3i];        
        if numel(p.rx.equalizer.StepSize_DDLMS) ~= 1 && numel(p.rx.equalizer.StepSize_DDLMS) ~= p.rx.equalizer.Iterations_DDLMS
            error('The DDLMS step size is incorrect. Please specify one step size for all iterations or a vector for the step size of each iteration.')
        end
        [~,ix]=sort(angle(p.Constellation));
        p.Constellation=p.Constellation(ix);
        p.rx.equalizer.Constellation= single(p.Constellation/sqrt(mean(abs(p.Constellation).^2)));          
    end
    
    % Storing ideal transmitted constellation
    try
        switch lower(p.modulation.format)
            case '16-qam'
                % 16QAM
%                 p.tx.ideal(k, :) = iq_modulator(p, ones(1, p.N_symb), [p.electric.levels(data_to_symnum(p.data(k).bits(1:2, :)));p.electric.levels(data_to_symnum(p.data(k).bits(3:4, :)))]);
%                 p.tx.ideal(k, :) = p.tx.ideal(k, :)./sqrt(mean(abs(p.tx.ideal(k, :)).^2)); % Normalizing
                dataSymbolsIn = bi2de(p.data(k).bits'); %Converting from bits to Symbols, note that p.data(k).bits is M*N where M is the number of bits per symbol.
                Sym_mod = qammod(dataSymbolsIn,16); %Modulating the symbols without oversampling.
                p.tx.ideal(k,:)=Sym_mod'./sqrt(mean(abs(Sym_mod').^2));
            case 'qpsk'
                % QPSK
                p.tx.ideal(k, :) = iq_modulator(p, ones(1, p.N_symb), [p.electric.levels(data_to_symnum(p.data(k).bits(1, :)));p.electric.levels(data_to_symnum(p.data(k).bits(2, :)))]);
                p.tx.ideal(k, :) = p.tx.ideal(k, :)./sqrt(mean(abs(p.tx.ideal(k, :)).^2)); % Normalizing
            case 'ook'
                % OOK
                p.tx.ideal(k, :) = iq_modulator(p, ones(1, p.N_symb), [p.electric.levels(data_to_symnum(p.data(k).bits(1, :)));0.5*ones(1,length(p.data(k).bits))]);
                p.tx.ideal(k, :) = p.tx.ideal(k, :)./sqrt(mean(abs(p.tx.ideal(k, :)).^2)); % Normalizing
                %p.tx.ideal(k,:)=p.data(k).bits(1, :);
            otherwise
                error('Modulation format %s not found', p.modulation.format)
        end
    catch
        error('Modulation format not specified (set p.modulation.format in get_parameters)')
    end
    
    % Set up the electric driving signals % TODO: Make polarization
    % compatible
    p.samp_per_symbol(k).wdm = p.samp_per_symb*(max(p.f_symb)/p.f_symb(k));
    switch lower(p.modulation.format)
        case '16-qam'
            %16QAM
            %p.electric.data = p.data(k).bits(1:2, 1:end/(max(p.f_symb)/p.f_symb(k)), 1); e_i = electric_signal(p,k);
            %p.electric.data = p.data(k).bits(3:4, 1:end/(max(p.f_symb)/p.f_symb(k)), 1); e_q = electric_signal(p,k);
            
            % Modulating data onto a new channel TODO: Make sure that all of the data is independent!!!
            %Filtering  the ideal signal with a RRC filter and oversampling.            
            span = 10;        % Filter span in symbols
            rolloff = 0.2;   % Roloff factor of filter
            txfilter = comm.RaisedCosineTransmitFilter('RolloffFactor',rolloff, 'FilterSpanInSymbols',span,'OutputSamplesPerSymbol',p.samp_per_symb); %Filter Object
            u_mod = txfilter(Sym_mod)';
        case 'qpsk'
            %QPSK
            p.electric.data = p.data(k).bits(1, 1:end/(max(p.f_symb)/p.f_symb(k)), 1); e_i = electric_signal(p,k);
            p.electric.data = p.data(k).bits(2, 1:end/(max(p.f_symb)/p.f_symb(k)), 1); e_q = electric_signal(p,k);
            % Modulating data onto a new channel TODO: Make sure that all of the data is independent!!!
            u_mod = iq_modulator(p, u_cw, [e_i; e_q]);
        case 'ook'
            %OOK
            p.electric.data = p.data(k).bits(1, 1:end/(max(p.f_symb)/p.f_symb(k)), 1); e_i = electric_signal(p,k);
            %p.electric.data = p.data(2, :, 1, k);
            e_q = 0.5*ones(1,length(e_i));
            % Modulating data onto a new channel TODO: Make sure that all of the data is independent!!!
            u_mod = iq_modulator(p, u_cw, [e_i; e_q]);            
    end
 
    %Seting the power to the desired output per channel.
    u_mod = set_power(u_mod, p.power_array(k));
    measure_power(u_mod, sprintf('after mod. (CH%d)', k));
    power_out = measure_power(u_mod);
    u_mod = u_mod.*exp(1j*p.w_index.sig(k)*p.t); % Shifting channel to new carrier frequency
    
    %Multipleing the WDM channels
    u0 = u0 + u_mod; % Adding new channel to the total output
    
    write_log(p, now, sprintf('Generated channel %d with %d symbols, Signal power after modulation, %.2f dBm, Center frequency, %.2f GHz (Relative to %.2f THz)', k, p.N_symb, power_out, p.w_index.sig(k)/(2*pi*1e9), p.const.c/p.const.lambda*1e-12));
end

% Convert to single precision; enables GPU acceleration, reduces disk-space
% of output files, but introduces an error that grows with the number of
% iterations in the split-step solver.
% Evaluate whether the extra precision is needed on a project-by-project
% basis.
% u0 = single(u0);

%p = rmfield(p, 'electric'); % Cleaning up useless fields

write_log(p, now, 'Signal Generation Complete');