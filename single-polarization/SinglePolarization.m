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
    
    %% Simulation parameters that can be set by users
    properties (Access=public)
        % User parallel in kmeans
        useParallel
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
            addParameter(p, 'linkArray', Link(), @(x) isa(x, 'Link'));
            addParameter(p, 'channelArray', Channel(), @(x) isa(x, 'Channel'));
            addParameter(p, 'useParallel', true, @islogical);
            
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
            obj.useParallel = p.Results.useParallel;
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
%                 fprintf(obj.logFid, 'Simulation Name, Simulation id\n');
                fprintf(obj.logFid, '%s, %d, %s\n', obj.simulationName, ...
                    obj.simulationId, datestr(now()));
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
            obj.N = ceil(4*obj.tmax*obj.fmax);
            
            obj.domega = 4*pi*obj.fmax/obj.N;
            obj.omega = (-obj.N/2:obj.N/2-1).'*obj.domega;
            
            obj.dt = 2*obj.tmax/obj.N;
            obj.t = (-obj.N/2:obj.N/2-1).'*obj.dt;
            
            %% Generate signals
            % Transmitted signal in time domain
            obj.txSignalTime = zeros(obj.N, 1);
            generateSignal(obj);
            % Transmitted signal in spectrum domain
            obj.txSignalSpectrum = ft(obj.txSignalTime, obj.domega);
            
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
        
        %% Simulate transmission and Receiver
        function simulate(obj)
            %% Transmission
            for n=1:obj.numberLink
                ssf(obj, n);
            end
            
            %% Receiver
            % Dispersion compensation
            dispersionCompensation(obj);
            for n=1:obj.numberChannel
                downConvert(obj, n);
                synchronize(obj, n);
                matchReceivedSignal(obj, n);
                findCloudCenter(obj, n);
                matchConstellation(obj, n);
                computeSER(obj, n);
                computeSNR(obj, n);
                computeEVM(obj, n);
            end
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
% spectrum. The factor of 10 leaves enough margin for the channel spectrum
channelSpectrumBound = abs([obj.channelArray.centerFrequency])+...
    symbolRate*10;

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
symbolTime = round(1e15./[obj.channelArray.symbolRate]);
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
% totalTime = totalTime*1.1;

% Compute actual number of symbols for every channel
actualNumberSymbol = ceil(totalTime./symbolTime);

% Compute tmax [s], convert back to second
tmax = totalTime/2/1e15;
end

function generateOOK(obj, channelIdx)
% Generate NRZ OOK without carrier suppress
channel = obj.channelArray(channelIdx);
assert (strcmp(channel.modulation, 'OOK'))

% Generate bit stream
channel.txBit = randi([0, 1], channel.actualNumberSymbol, channel.bitPerSymbol);
% The first and last ceil(channel.symbolInFir/2) symbols are affected by
% FIR convolution effects. Set them to 0 to eliminate this effect.
halfFirSymbolLength = ceil(channel.symbolInFir/2);
channel.txBit(1:halfFirSymbolLength) = 0;
channel.txBit(end-halfFirSymbolLength+1:end) = 0;

% Symbol stream
channel.txSymbol = channel.txBit;
channel.txTime = repmat(channel.txBit, 1, channel.actualSamplePerSymbol).';
channel.txTime = channel.txTime(:);
dataTimeLength = obj.N;

% Generate Gauss FIR
channel.fir = gaussdesign(channel.firFactor, channel.symbolInFir, channel.actualSamplePerSymbol);
% Pass txTime through FIR
channel.txTime = upfirdn(channel.txTime, channel.fir);

% Shift signal randomly within one symbol time
channel.shiftNumberSample = randi([0, channel.actualSamplePerSymbol-1]);
channel.txTime = circshift(channel.txTime, channel.shiftNumberSample);

% FIR delay in number of samples
firOverhead = ceil((length(channel.fir)-1)/2);
% Remove head and tail added by FIR and convolution inside upfirdn
s = (firOverhead+1):(firOverhead+dataTimeLength);
channel.txTime = channel.txTime(s);
end

function generate16QAM(obj, channelIdx)
% Generate 16QAM bits
channel = obj.channelArray(channelIdx);
assert (strcmp(channel.modulation, '16QAM'))

channel.txBit = randi([0, 1], channel.actualNumberSymbol, channel.bitPerSymbol);
% The first and last ceil(channel.symbolInFir/2) symbols are affected by
% FIR convolution effects. Set them to 0 to eliminate this effect.
halfFirSymbolLength = ceil(channel.symbolInFir/2);
channel.txBit(1:halfFirSymbolLength, :) = 0;
channel.txBit(end-halfFirSymbolLength+1:end, :) = 0;

% Symbols
channel.txSymbol = bi2de(channel.txBit);
channel.txSymbol = qammod(channel.txSymbol, channel.constellationSize);
channel.txSymbol = channel.txSymbol/sqrt(mean(abs(channel.txSymbol).^2)); %

% Square root raised cosine FIR
channel.fir = rcosdesign(channel.firFactor, channel.symbolInFir, channel.actualSamplePerSymbol, 'sqrt');
% Length of txTime
dataTimeLength = obj.N;
% Pass through FIR with upsample
channel.txTime = upfirdn(channel.txSymbol, channel.fir, channel.actualSamplePerSymbol);

% Shift signal randomly within one symbol time
channel.shiftNumberSample = randi([0, channel.actualSamplePerSymbol-1]);
channel.txTime = circshift(channel.txTime, channel.shiftNumberSample);

% Remove head and tail of txTime
% This is the head to be removed, because the output length of FIR is
% ceil(((length(xin)-1)*p+length(h))/q) for yout = upfirdn(xin,h,p,q)
% See Matlab reference of upfirdn and conv
firOverhead = ceil((length(channel.fir)-channel.actualSamplePerSymbol)/2);
s = (firOverhead+1):(firOverhead+dataTimeLength);
channel.txTime = channel.txTime(s);
end

function generateSignal(obj)
% Generate signal according to modulation format

for n=1:obj.numberChannel
    channel = obj.channelArray(n);
    if strcmp(channel.modulation, 'OOK')
        generateOOK(obj, n);
    elseif strcmp(channel.modulation, '16QAM')
        generate16QAM(obj, n);
    end
    
    % Assign power to the signal
    powerNormalized = norm(channel.txTime)^2/obj.N;
    channel.txTime = channel.txTime*sqrt(channel.powerW/powerNormalized);
    
    % Move in spectrum to its center frequency
    channel.txTime = channel.txTime.*...
        exp(-1i*2*pi*channel.centerFrequency.*obj.t);
    
    % Add channel signal to the total transmitted signal
    obj.txSignalTime = obj.txSignalTime+channel.txTime;
end
end

function xf = ft(xt, domega)
% Fourier transform
xf = fftshift(ifft(fftshift(xt)))/sqrt(domega/(2*pi));
end

function xt = ift(xf, domega)
% Inverse Fourier transform
xt = fftshift(fft(fftshift(xf)))*sqrt(domega/(2*pi));
end

function ssf(obj, linkIdx)
% Split step Fourier

link = obj.linkArray(linkIdx);

%% Dispersive and nonlinear phase factors
dispersion2 = exp(0.5*1i*link.beta2*obj.omega.^2*link.dz);
dispersion3 = exp(1i/6*link.beta3*obj.omega.^3*link.dz);
dispersion = dispersion2.*dispersion3;
hhz = 1i*link.gamma*link.dzEff;

%% Main loop
% scheme: 1/2N -> D -> 1/2N; first half step nonlinear
uu = obj.currentSignalTime; % time domain signal
temp = uu.*exp(abs(uu).^2.*hhz/2); % note hhz/2
for n=1:link.numberSteps
    % dispersion
    uu = ift(ft(temp, obj.domega).*dispersion, obj.domega);
    % nonlinearity
    temp = uu.*exp(abs(uu).^2.*hhz).*exp(-0.5*link.alphaLinear*link.dz);
end
uu = temp.*exp(-abs(uu).^2.*hhz/2); % Final field

%% EDFA
uu = uu*exp(0.5*link.alphaLinear*link.spanLength);

% Power spectral density (PSD) of ASE noise, G=(exp(alpha*L)-1)*h*nu*nsp
noisePower = (exp(link.alphaLinear*link.spanLength))*link.h*link.nu*link.nsp;
% Convert PSD to signal power in the time domain
noisePower = noisePower*obj.N*obj.domega/(2*pi);
% add noise
uu = uu + sqrt(0.5*noisePower)*(randn(size(uu, 1), 1)+1i*randn(size(uu, 1), 1));

%% DCF
% DCF is an ideal FBG
% Note the signs of beta2 and beta3 are changed to compensate for
% dispersion
fbg = exp(-0.5*1i*link.beta2*obj.omega.^2*link.DCFLength).*...
    exp(-1i/6*link.beta3*obj.omega.^3*link.DCFLength);

temp = ft(uu, obj.domega).*fbg;
uu = ift(temp, obj.domega);

%% Output signal in time and spectrum domains
obj.currentSignalTime = uu;
obj.currentSignalSpectrum = temp;

end

function dispersionCompensation(obj)
% Compute residual dispersion eqivalent length

fbg = ones(size(obj.omega));
for n=1:obj.numberLink
    link = obj.linkArray(n);
    residualDispersionLength = link.spanLength-link.DCFLength;
    fbg = fbg.*exp(-0.5*1i*link.beta2*obj.omega.^2*residualDispersionLength).*...
        exp(-1i/6*link.beta3*obj.omega.^3*residualDispersionLength);
end

obj.currentSignalSpectrum = obj.currentSignalSpectrum.*fbg;
obj.currentSignalTime = ift(obj.currentSignalSpectrum, obj.domega);

end

function downConvert(obj, channelIdx)
% DSP at receiver
signal = obj.currentSignalTime;
channel = obj.channelArray(channelIdx);

% Down convert
signal = signal.*exp(1i*2*pi*channel.centerFrequency.*obj.t);
signal = upfirdn(signal, channel.fir, 1, 1);

% Compute overhead
overhead = ceil((channel.fir-1)/2);
signal = signal((overhead+1):(overhead+obj.N));

channel.rxTime = signal;
end

function synchronize(obj, channelIdx)
% Find the optimal place to sample the signal

channel = obj.channelArray(channelIdx);
% Copy received signal
rxSignal = channel.rxTime;
% Real and imaginary parts of the signal
signal = zeros(size(rxSignal, 1), 2);
signal(:, 1) = real(rxSignal);
signal(:, 2) = imag(rxSignal);

% Sum of distances from points to their centers
sumDistance = zeros(channel.actualSamplePerSymbol, 1);
% Sum of signal to noise ratios of all point clouds
q = zeros(channel.actualSamplePerSymbol, 1);
% Another measurement
q2 = zeros(channel.actualSamplePerSymbol, 1);

% Exhaustively try each possible sampling offset and calculate the two
% measurements
for nn=1:channel.actualSamplePerSymbol
    % Downsample with the given offset
    tmpSignal = downsample(signal, channel.actualSamplePerSymbol, nn-1);
    % Instruct kmeans to use parallel threads
    opts = statset('UseParallel', obj.useParallel);
    [~, tmpCenters, tmpDistances] = ...
        kmeans(tmpSignal, channel.constellationSize, ...
        'Display', 'off', 'maxiter', 1000, ...
        'Replicates', 4, 'Options', opts);
    
    sumDistance(nn) = sum(tmpDistances);
    
    tmpSignalPower = tmpCenters(:, 1)+1i*tmpCenters(:, 2);
    tmpSignalPower0 = abs(tmpSignalPower-tmpSignalPower.');
    
    tmpSignalPower = tmpSignalPower0;
    tmpSignalPower(tmpSignalPower==0) = inf;
    tmpSignalPower = min(tmpSignalPower(:));
    q(nn) = tmpSignalPower/sumDistance(nn);
    
    tmpq2 = tmpSignalPower0./(tmpDistances+tmpDistances.');
    q2(nn) = mean(tmpq2(:));
end
[~, sumDistanceIdx] = min(sumDistance);
[~, qIdx] = max(q);
[~, q2Idx] = max(q2);
channel.rxOptimalOffset = q2Idx;

% Down sample received signal at the optimal offset
signal = downsample(signal, channel.actualSamplePerSymbol, ...
    channel.rxOptimalOffset-1);
channel.rxSymbol = signal(:, 1)+1i*signal(:, 2);
end

function matchReceivedSignal(obj, channelIdx)
% Match the received signal with the transmitted signal, so that later on
% we can calculate SER, SNR, EVM, and so on.

channel = obj.channelArray(channelIdx);
% Transmitted symbols
txSymbol = channel.txSymbol;
% Received symbols
rxSymbol = channel.rxSymbol;
% Correlation between received and transmitted symbols
[correlation, lag] = xcorr(rxSymbol, txSymbol);
% Find the offset of the maximum correlation
[~, maxIdx] = max(correlation);
maxLag = lag(maxIdx);

% if maxLag<0, rxSymbol remove abs(maxLag) tail elements and txSymbol
% remove abs(maxLag) head elements
% if maxLag>0, rxSymbol remove abs(maxLag) head elements and txSymbol
% remove abs(maxLag) tail elements
% if maxLag==0, keep them the same
if maxLag<0
    rxSymbol = rxSymbol(1:end+maxLag);
    txSymbol= txSymbol(1-maxLag:end);
elseif maxLag>0
    rxSymbol = rxSymbol(1+maxLag:end);
    txSymbol = txSymbol(1:end-maxLag);
end
channel.rxSymbolMatched = rxSymbol;
channel.txSymbolMatched = txSymbol;
end

function findCloudCenter(obj, channelIdx)
% Find centers of point clouds
opts = statset('UseParallel', obj.useParallel);
symbol = obj.channelArray(channelIdx).rxSymbolMatched;
symbol = [real(symbol), imag(symbol)];
[~, center] = kmeans(symbol, ...
    obj.channelArray(channelIdx).constellationSize, ...
    'Display', 'off', 'maxiter', 1000, ...
    'Replicates', 16, 'Options', opts);
obj.channelArray(channelIdx).rxCloudCenter = center(:, 1)+1i*center(:, 2);
end

function matchConstellation(obj, channelIdx)
% match received and transmitted
tx = obj.channelArray(channelIdx).txSymbolMatched;
rx = obj.channelArray(channelIdx).rxSymbolMatched;

% Objective function to minimize for each pair of tx and rx points
% objfcn = @(x) x(1)+x(2)*rx - tx;
% x0 = (1+1i)*[1;1]; % arbitrary initial guess
rxV = zeros(length(rx), 2);
rxV(:, 1) = real(rx);
rxV(:, 2) = imag(rx);
txV = zeros(length(tx), 2);
txV(:, 1) = real(tx);
txV(:, 2) = imag(tx);

objfcn = @(x) sum((rxV(:, 1)*x(1) - rxV(:, 2)*x(2) - txV(:, 1)).^2 + ...
    (rxV(:, 2)*x(1) + rxV(:, 1)*x(2) - txV(:, 2)).^2);
x0 = [1; 1];
% Optimize nonlinear least-squares problem
opts = optimoptions(@fminunc,'Display', 'off', ...
'UseParallel', obj.useParallel, 'MaxFunctionEvaluations', 1e4);

% Solution 
[xSol,fval,exitflag,output] = fminunc(objfcn,x0,opts);
xSol = xSol(1)+1i*xSol(2);
% 
% obj.channelArray(channelIdx).rxSymbolRotated = xSol(1) + xSol(2)*rx;
obj.channelArray(channelIdx).rxSymbolRotated = xSol*rx;
end

function computeEVM(obj, channelIdx)

powerVector = abs(obj.channelArray(channelIdx).txSymbolMatched).^2;
errorVector = abs(obj.channelArray(channelIdx).txSymbolMatched-...
    obj.channelArray(channelIdx).rxSymbolRotated).^2;
obj.channelArray(channelIdx).EVM = ...
    sqrt(mean(errorVector)/mean(powerVector));
end

function computeSER(obj, channelIdx)

tx = obj.channelArray(channelIdx).txSymbolMatched;
% Constellation points of transmitted symbols
txUnique = unique(tx);
% Constellation points of received symbols
rx = obj.channelArray(channelIdx).rxSymbolRotated;
% Distances between received symbol to all transmitted symbols
dist = abs(rx-txUnique.');
% For each received symbol, calculate its distance to the closes
% transmitted symbol
[minDist, ~] = min(dist, [], 2);
% For each received symbol, calculate its distance to its corresponding
% transmitted symbol
txDist = abs(rx-tx);
% Decode a symbol to its nearest transmitted symbol, and calculate the
% corresponding SER
obj.channelArray(channelIdx).SER = sum(minDist<txDist)/length(txDist);
end

function computeSNR(obj, channelIdx)
rx = obj.channelArray(channelIdx).rxSymbolRotated;
tx = obj.channelArray(channelIdx).txSymbolMatched;
txUnique = unique(tx);

powerSignal = 0;
powerNoise = 0;
% Iterate over all unique transmitted symbols
for n=1:length(txUnique)
    tmpSignal = rx(tx==txUnique(n));
    tmpRatio = length(tmpSignal)/length(rx);
    powerSignal = powerSignal + abs(mean(tmpSignal))^2*tmpRatio;
    powerNoise = powerNoise + abs(std(tmpSignal))^2*tmpRatio;
end
obj.channelArray(channelIdx).SNR = powerSignal/powerNoise;
obj.channelArray(channelIdx).SNRdB = ...
    10*log10(obj.channelArray(channelIdx).SNR);
end