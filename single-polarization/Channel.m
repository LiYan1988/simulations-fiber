classdef Channel
    %Optical channel
    %   Detailed explanation goes here
    
    properties
        modulation % {'OOK', '16QAM'}
        symbolRate % [Hz]
        centerFrequency % [Hz]
        minSamplePerSymbol % minimum sample per symbol, can be higher
        powerdB % power [dB]
        powerW % power [W]
        
        % Filter parameters
        filterFactor % bt of gaussdesign or beta of rcosdesign
        symbolInFilter % number of symbols in FIR
    end
    
    properties (Dependent)
        constellationSize % number of points in the constellation diagram
    end
    
    methods
        function obj = Channel(varargin)
            %Construct an instance of Channel
            %   Inputs are name-value pairs, default value is for OOK
            
            %% Parse input
            p = inputParser;
            
            validModulation = {'OOK', '16QAM'};
            checkModulation = @(x)any(validatestring(x, validModulation));
            addParameter(p, 'modulation', 'OOK', checkModulation);
            addParameter(p, 'symbolRate', 10e9, @isnumeric);
            addParameter(p, 'centerFrequency', 0e9, @isnumeric);
            addParameter(p, 'minSamplePerSymbol', 4, @isnumeric);
            addParameter(p, 'powerdB', 0, @isnumeric);
            
            % Filter parameters
            addParameter(p, 'filterFactor', 0.7, @isnumeric);
            addParameter(p, 'symbolInFilter', 10, @isnumeric);
            
            % parse the inputs
            parse(p, varargin{:});
            
            %% Set parameters
            obj.modulation = p.Results.modulation;
            obj.symbolRate = p.Results.symbolRate;
            obj.centerFrequency = p.Results.centerFrequency;
            obj.minSamplePerSymbol = p.Results.minSamplePerSymbol;
            obj.powerdB = p.Results.powerdB;
            obj.powerW = 10.^(obj.powerdB/10)/1e3;
            
            obj.filterFactor = p.Results.filterFactor;
            obj.symbolInFilter = p.Results.symbolInFilter;
        end
        
        function constellationSize = get.constellationSize(obj)
            if strcmp(obj.modulation, 'OOK')
                constellationSize = 2;
            elseif strcmp(obj.modulation, '16QAM')
                constellationSize = 16;
            end
        end
        
        
    end
end

