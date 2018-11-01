classdef SinglePolarization
    %Single polarization simulation of hybrid channels
    %   Detailed explanation goes here
    
    %% Constant simulation data
    properties(Constant)
        % [m/s], Speed of light
        lightSpeed = 2.99792458e8;
        % [m], reference wavelength
        wavelength = 1550e-9;
        % [J*s] or [W*Hz^-2], Plank's constant
        h = 6.626e-34;
        % [Hz], frequency
        nu = 1.934144890322581e+14;
    end
    
    %% Regular simulation data
    properties
        % Simulation Name
        simulationName
        % Simulation id
        simulationId
        % Random seed
        randomSeed = 0;
        % Result folder name
        resultFolder
        % Log file name
        logFile
        % Log file fid
        logFid
        
        % Links
        linkArray
        
        % Channels
        channelArray
    end
    
    properties (Dependent)
        numberLink
        numberChannel
    end
    
    methods
        function obj = SinglePolarization(varargin)
            %Construct an instance of this class
            %   Inputs are name-value pairs
            
            %% Parse input
            p = inputParser;
            
            addParameter(p, 'simulationName', 'singlePolarization', @ischar);
            addParameter(p, 'simulationId', 0, @isnumeric);
            addParameter(p, 'randomSeed', 0, @isnumeric);
            addParameter(p, 'resultFolder', 'results/', @ischar);
            addParameter(p, 'logFile', '', @ischar);
            addParameter(p, 'linkArray', Link(), @(x) isa(x, 'Link'));
            addParameter(p, 'channelArray', Channel(), @(x) isa(x, 'Channel'));
            
            % Parse inputs
            parse(p, varargin{:});
            
            %% Set parameters
            obj.simulationName = p.Results.simulationName;
            obj.simulationId = p.Results.simulationId;
            obj.randomSeed = p.Results.randomSeed;
            obj.resultFolder = p.Results.resultFolder;
            obj.logFile = sprintf('%s_%d.csv', obj.simulationName, obj.simulationId);
            obj.linkArray = p.Results.linkArray;
            obj.channelArray = p.Results.channelArray;
            
            %% Setup result folder and log file
            % Result folder
            if ~exist(obj.resultFolder, 'dir')
                mkdir(obj.resultFolder)
            end
            % Log file
            obj.logFile = fullfile(obj.resultFolder, obj.logFile);
            
            % Open log file
            obj.logFid = fopen(obj.logFile, 'a');
            if obj.logFid == -1
                warning('Could not open log file');
            else
                fprintf(obj.logFid, 'Simulation Name, Simulation id\n');
                fprintf(obj.logFid, '%s, %d\n', obj.simulationName, obj.simulationId);
            end
        end
        
        function numberLink = get.numberLink(obj)
            numberLink = length(obj.linkArray);
        end
        
        function numberChannel = get.numberChannel(obj)
            numberChannel = length(obj.channelArray);
        end
        
        function obj = sortChannelArray(obj)
            % Sort channelArray, reorder channelArray according to their
            % center frequencies
            freq = [obj.channelArray.centerFrequency];
            [~, idx] = sort(freq);
            obj.channelArray = obj.channelArray(idx);
        end
        
    end
end

%% Helper Functions for Class Methods
