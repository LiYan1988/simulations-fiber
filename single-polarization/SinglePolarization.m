classdef SinglePolarization < handle
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
    properties (GetAccess=public, SetAccess=private)
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
        
        % Max spectrum bandwidth [Hz], the total spectrum ranges from -fmax
        % to fmax
        fmax
        % Max time span [s], the total time ranges from -tmax to tmax
        tmax
        % Number of samples in simulation
        N
        % delta omega [rad/s]
        domega
        % omega axis [rad/s], a column vector
        omega
        % delta time [s]
        dt
        % time axis [s], a column vector
        t
        
        % Transmitted signal in time domain
        txSignalTime
        % Transmitted signal in spectrum domain
        txSignalSpectrum
        % Current signal in time domain
        currentSignalTime
        % Current Signal in spectrum domain
        currentSignalSpectrum
    end
    
    %% Dependent parameters
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
            obj.linkArray = copy(p.Results.linkArray);
            obj.channelArray = copy(p.Results.channelArray);
            
            %% Set random seed for Matlab
            rng(obj.randomSeed)
            
            %% Sort channel array
            freq = [obj.channelArray.centerFrequency];
            [~, idx] = sort(freq);
            obj.channelArray = obj.channelArray(idx);
            
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
            
            %% Calculate fmax
            % The spectrum axis ranges from -fmax to fmax
            [obj.fmax, actualSamplePerSymbol] = computeTotalSpectrum(obj);
            for n=1:obj.numberChannel
                obj.channelArray(n).actualSamplePerSymbol = ...
                    actualSamplePerSymbol(n);
            end
            
            %% Calculate tmax
            % The time axis ranges from -tmax to tmax
            [obj.tmax, actualNumberSymbol] = computeTotalTime(obj);
            for n=1:obj.numberChannel
                obj.channelArray(n).actualNumberSymbol = ...
                    actualNumberSymbol(n);
            end
            
            %% Create spectrum and time axis
            obj.N = 4*obj.tmax*obj.fmax;
            
            obj.domega = 4*pi*obj.fmax/obj.N;
            obj.omega = (-obj.N/2:obj.N/2-1).'*obj.domega;
            
            obj.dt = 2*obj.tmax/obj.N;
            obj.t = (-obj.N/2:obj.N/2-1).'*obj.dt;
            
            %% Generate signals
            % Transmitted signal in time domain
            obj.txSignalTime = zeros(obj.N, 1);
            for n=1:obj.numberChannel
                generateSignal(obj.channelArray(n), obj);
            end
            % Transmitted signal in spectrum domain
            obj.txSignalSpectrum = ft(obj.txSignalTime, obj.dt);
            
            % Current signal in time and spectrum domains
            obj.currentSignalTime = obj.txSignalTime;
            obj.currentSignalSpectrum = obj.txSignalSpectrum;
        end
        
        function numberLink = get.numberLink(obj)
            numberLink = length(obj.linkArray);
        end
        
        function numberChannel = get.numberChannel(obj)
            numberChannel = length(obj.channelArray);
        end
    end
end

%% Helper Functions for Class Methods
function [fmax, actualSamplePerSymbol] = computeTotalSpectrum(obj)
% Compute the total spectrum bandwidth for simulation
% Inputs:
%   obj: SinglePolarization object
% Outputs:
%   fmax: the single side bandwidth of the total spectrum, the overall
%       total bandwidth is 2*fmax
%   actualSamplePerSymbol: actual sample per symbol when fmax is
%       computed

symbolRate = [obj.channelArray.symbolRate];
minSamplePerSymbol = [obj.channelArray.minSamplePerSymbol];
% Bounds of channel spectrum, which should be covered by the half total
% spectrum
channelSpectrumBound = abs([obj.channelArray.centerFrequency])+...
    symbolRate/2;

% Compute the least common multiple of all symbol rates
totalSpectrum = 1;
for n=1:length(symbolRate)
    totalSpectrum = lcm(totalSpectrum, symbolRate(n));
end

% Check if the sample per symbol is greater than the minimum
% for all the channels. If not, multiply with 2 until it satisfies the
% requirement.
for n=1:length(minSamplePerSymbol)
    while totalSpectrum/symbolRate(n)<minSamplePerSymbol(n)
        totalSpectrum = 2*totalSpectrum;
    end
end

% Check that all the channels are within the total spectrum
for n=1:length(channelSpectrumBound)
    while totalSpectrum/2<channelSpectrumBound(n)
        totalSpectrum = 2*totalSpectrum;
    end
end

% Compute actual sample per symbol
actualSamplePerSymbol = totalSpectrum./symbolRate;

% fmax
fmax = totalSpectrum/2;
end

function [tmax, actualNumberSymbol] = computeTotalTime(obj)
% Compute the total time span for simulation
% Inputs:
%   obj: SinglePolarization object
% Outputs:
%   tmax: the single side time span of the simulation, the overall
%       total time span is 2*tmax
%   actualNumberSymbol: actual number of symbols when tmax is computed

% Symbol time [fs], convert to fs so that all symbol times are integer
symbolTime = 1e15./[obj.channelArray.symbolRate];
minNumberSymbol = [obj.channelArray.minNumberSymbol];

% Compute the least common multiple of all symbol times
totalTime = 1;
for n=1:length(symbolTime)
    totalTime = lcm(totalTime, symbolTime(n));
end

% Check if the total symbol time is greater than the minimum
% for all the channels. If not, multiply with 2 untile it satisfies all
% the requirements.
for n=1:length(minNumberSymbol)
    while totalTime/symbolTime(n)<minNumberSymbol(n)
        totalTime = 2*totalTime;
    end
end

% Compute actual number of symbols for every channel
actualNumberSymbol = totalTime./symbolTime;

% Compute tmax [s], convert back to second
tmax = totalTime/2/1e15;
end

function generateOOK(channel)
% Generate NRZ OOK without carrier suppress

% Generate bit stream
channel.dataBit = randi([0, 1], channel.actualNumberSymbol, channel.bitPerSymbol);
% The first and last ceil(channel.symbolInFir/2) symbols are affected by
% FIR convolution effects. Set them to 0 to eliminate this effect.
halfFirSymbolLength = ceil(channel.symbolInFir/2);
channel.dataBit(1:halfFirSymbolLength) = 0;
channel.dataBit(end-halfFirSymbolLength+1:end) = 0;

% Symbol stream
channel.dataSymbol = channel.dataBit;
channel.dataTime = repmat(channel.dataBit, 1, channel.actualSamplePerSymbol).';
channel.dataTime = channel.dataTime(:);
dataTimeLength = length(channel.dataTime);

% Generate Gauss FIR
channel.fir = gaussdesign(channel.firFactor, channel.symbolInFir, channel.actualSamplePerSymbol);
% Pass dataTime through FIR
channel.dataTime = upfirdn(channel.dataTime, channel.fir);

% Shift signal randomly within one symbol time
channel.shiftNumberSample = randi([0, channel.actualSamplePerSymbol-1]);
% channel.dataTime = [zeros(channel.shiftNumberSample, 1); ...
%     channel.dataTime(1:end-channel.shiftNumberSample)];
channel.dataTime = circshift(channel.dataTime, channel.shiftNumberSample);

% FIR delay in number of samples
firOverhead = ceil((length(channel.fir)-1)/2);
% Remove head and tail added by FIR and convolution inside upfirdn
s = (firOverhead+1):(firOverhead+dataTimeLength);
channel.dataTime = channel.dataTime(s);
end

function generate16QAM(channel)
% Generate 16QAM bits
channel.dataBit = randi([0, 1], channel.actualNumberSymbol, channel.bitPerSymbol);
% The first and last ceil(channel.symbolInFir/2) symbols are affected by
% FIR convolution effects. Set them to 0 to eliminate this effect.
halfFirSymbolLength = ceil(channel.symbolInFir/2);
channel.dataBit(1:halfFirSymbolLength, :) = 0;
channel.dataBit(end-halfFirSymbolLength+1:end, :) = 0;

% Symbols
channel.dataSymbol = bi2de(channel.dataBit);
channel.dataSymbol = qammod(channel.dataSymbol, channel.constellationSize);
channel.dataSymbol = channel.dataSymbol/sqrt(mean(abs(channel.dataSymbol).^2)); %

% Square root raised cosine FIR
channel.fir = rcosdesign(channel.firFactor, channel.symbolInFir, channel.actualSamplePerSymbol, 'sqrt');
% Length of dataTime
dataTimeLength = channel.actualSamplePerSymbol*channel.actualNumberSymbol;
% Pass through FIR with upsample
channel.dataTime = upfirdn(channel.dataSymbol, channel.fir, channel.actualSamplePerSymbol);

% Shift signal randomly within one symbol time
channel.shiftNumberSample = randi([0, channel.actualSamplePerSymbol-1]);
% channel.dataTime = [zeros(channel.shiftNumberSample, 1); ...
%     channel.dataTime(1:end-channel.shiftNumberSample)];
channel.dataTime = circshift(channel.dataTime, channel.shiftNumberSample);

% Remove head and tail of dataTime
% This is the head to be removed, because the output length of FIR is
% ceil(((length(xin)-1)*p+length(h))/q) for yout = upfirdn(xin,h,p,q)
% See Matlab reference of upfirdn and conv
firOverhead = ceil((length(channel.fir)-channel.actualSamplePerSymbol)/2);
s = (firOverhead+1):(firOverhead+dataTimeLength);
channel.dataTime = channel.dataTime(s);
end

function generateSignal(channel, obj)
% Generate signal according to modulation format
if strcmp(channel.modulation, 'OOK')
    generateOOK(channel);
elseif strcmp(channel.modulation, '16QAM')
    generate16QAM(channel);
end

% Assign power to the signal
powerNormalized = norm(channel.dataTime)^2/obj.N;
channel.dataTime = channel.dataTime*sqrt(channel.powerW/powerNormalized);

% Move in spectrum to its center frequency
channel.dataTime = channel.dataTime.*...
    exp(-1i*2*pi*channel.centerFrequency.*obj.t);

% Add channel signal to the total transmitted signal
obj.txSignalTime = obj.txSignalTime+channel.dataTime;
end

function xf = ft(xt, df)
% Fourier transform
xf = fftshift(ifft(fftshift(xt)))/sqrt(df/(2*pi));
end

function xt = ift(xf, df)
% Inverse Fourier transform
xt = fftshift(fft(fftshift(xf)))*sqrt(df/(2*pi));
end