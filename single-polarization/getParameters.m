function [param, inputs] = getParameters(varargin)
% Set up parameters for simulation.
%
% INPUTS:
%     linkConfig: a cell containing parameters for each link, it has higher
%     priority than numberLink. It is a cell of the following name-value
%     pairs, use param.basicLink as the default value if not stated in the cell.
%                 --------------------------------------------
%                 |          Name: Example values            |
%                 |    spanLength: 82000                     |
%                 |             D: 1.700000000000000e-05     |
%                 |         beta2: -2.168261939141489e-26    |
%                 |             S: -21.935483870967740       |
%                 |         beta3: 0                         |
%                 |     alphaIndB: 2.000000000000000e-04     |
%                 |   alphaLinear: 4.605170185988093e-05     |
%                 |         gamma: 0.001400000000000         |
%                 |   numberSteps: 100                       |
%                 |          NFdB: 5.563000000000000         |
%                 |           nsp: 1.799989635216783         |
%                 |     DCFLength: 82000                     |
%                 --------------------------------------------
% 
% OUTPUTS:
%   p: The parameter struct
%
% param = getParameters(
% 'randomSeed', 1, 
% 'resultsFolder', 'testResults', 
% 'logFileName', 'testLog', 
% 'numberSpan', 10, 
% 'spanLength', 82e3, 
% 'D', 17e-6, 
% 'DCFLength', 82e3
% );
% 

%% Parse inputs
% Create an InputParser object
p = inputParser;

% Input check functions
checkPositiveScalar = @(x, variableName) assert(isnumeric(x) && ...
    isscalar(x) && (x>0), ...
    sprintf('%s should be a positive number.', variableName));

checkScalar = @(x, variableName) assert(isnumeric(x) && isscalar(x), ...
    sprintf('%s should be a number', variableName));

checkCell = @(x, variableName) assert(iscell(x),...
    sprintf('%s should be a cell.', variableName));


% Add input schemes
% Random seed
defaultRandomSeed = 0;
checkRandomSeed = @(x) checkPositiveScalar(x, 'Random Seed');
addParameter(p, 'randomSeed', defaultRandomSeed, checkRandomSeed)

% Results folder
defaultResultsFolder = 'results/';
addParameter(p, 'resultsFolder', defaultResultsFolder, @ischar)

% Log file name
defaultLogFileName = '';
addParameter(p, 'logFileName', defaultLogFileName, @ischar);

% Nubmer of spans
defaultNumberSpan = 1;
checkNumberSpan = @(x) checkPositiveScalar(x, 'Number of spans');
addParameter(p, 'numberSpan', defaultNumberSpan, checkNumberSpan)

% Span length
defaultSpanLength = 82e3; % [m]
checkSpanLength = @(x) checkPositiveScalar(x, 'Span length');
addParameter(p, 'spanLength', defaultSpanLength, checkSpanLength);

% DCF length
defaultDCFLength = 82e3; % [m]
checkDCFLength = @(x) checkPositiveScalar(x, 'DCF length');
addParameter(p, 'DCFLength', defaultDCFLength, checkDCFLength);

% Dispersion parameter [s/m^2], correspond to beta2 (GVD) [s^2/m]
defaultD = 17e-6;
checkD = @(x) checkScalar(x, 'D');
addParameter(p, 'D', defaultD, checkD);

% S: Slop of D [s/m^3], correspond to beta3 (TOD) [s^3/m]
defaultS = -21.935483870967740; % default value corresponds to beta3=0
checkS = @(x) checkScalar(x, 'S');
addParameter(p, 'S', defaultS, checkS);

% Attenuation [dB/m], positive number
defaultAlphaIndB = 0.2e-3;
checkAlphaIndB = @(x) checkPositiveScalar(x, 'Alpha in dB');
addParameter(p, 'alphaIndB', defaultAlphaIndB, checkAlphaIndB);

% Kerr nonlinearity [1/W/m]
defaultGamma = 1.4e-3;
checkGamma = @(x) checkPositiveScalar(x, 'Gamma');
addParameter(p, 'gamma', defaultGamma, checkGamma);

% Number of steps
defaultNumberSteps = 100;
checkNumberSteps = @(x) checkPositiveScalar(x, 'Number of SSF steps');
addParameter(p, 'numberSteps', defaultNumberSteps, checkNumberSteps);

% Noise figure, NF in dB
defaultNFdB = 5.563;
checkNFdB = @(x) checkPositiveScalar(x, 'Noise figure in dB');
addParameter(p, 'NFdB', defaultNFdB, checkNFdB);

% Number of links
defaultnumberLink = 1;
checknumberLink = @(x) checkPositiveScalar(x, 'Number of links');
addParameter(p, 'numberLink', defaultnumberLink, checknumberLink);

% Link configuration
defaultChannelConfig = cell(1);
checkLinkConfig = @(x) checkCell(x, 'LinkConfig');
addParameter(p, 'linkConfig', defaultChannelConfig, checkLinkConfig);

% Modulation
defaultModulation = 'OOK';
validModulation = {'OOK', '16QAM'};
addParameter(p, 'modulation', defaultModulation, ...
    @(x)any(validatestring(x, validModulation)));

% Symbol rate [Hz]
defaultSymbolRate = 10e9;
checkSymbolRate = @(x) checkPositiveScalar(x, 'Symbol rate');
addParameter(p, 'symbolRate', defaultSymbolRate, checkSymbolRate);

% Center frequency [Hz]
defaultCenterFrequency = 0;
checkCenterFrequency = @(x) checkScalar(x, 'Center frequency');
addParameter(p, 'centerFrequency', ...
    defaultCenterFrequency, checkCenterFrequency);

% Filter parameter for Gaussian or square-root raised cosine FIR
defaultFilterFactor = 0.7;
checkFilterFactor = @(x) checkPositiveScalar(x, 'Channel filter factor');
addParameter(p, 'filterFactor', defaultFilterFactor, checkFilterFactor);

% Min sample per symbol
defaultMinSamplePerSymbol = 4;
checkMinSamplePerSymbol = @(x) checkPositiveScalar(x, ...
    'Minimum sample per symbol');
addParameter(p, 'minSamplePerSymbol', defaultMinSamplePerSymbol, ...
    checkMinSamplePerSymbol);

% Number of channels
defaultNumberChannel = 1;
checkNumberChannel = @(x) checkPositiveScalar(x, 'Number of channels');
addParameter(p, 'numberChannel', defaultNumberChannel, checkNumberChannel);

% Channel configuration
defaultChannelConfig = cell(1);
checkChannelConfig = @(x) checkCell(x, 'ChannelConfig');
addParameter(p, 'channelConfig', defaultChannelConfig, checkChannelConfig);

% Parse the inputs
parse(p, varargin{:})

%% Flags
param.flag.dualPol = 0; % Flag for dual-polarization

%% Constants
% [m/s], Speed of light
param.constant.lightSpeed = 2.99792458e8;
% [m], reference wavelength
param.constant.wavelength = 1550e-9;
% [J*s] or [W*Hz^-2], Plank's constant
param.constant.h = 6.626e-34;
% [Hz], frequency
param.constant.nu = param.constant.lightSpeed/param.constant.wavelength;

%% Random seed
% Random seed
param = setfield(param, 'random', 'seed', p.Results.randomSeed);

%% Log file for results
% In the same directory, create a folder and save logs into it
param = setfield(param, 'result', 'folder', p.Results.resultsFolder);
if ~exist(param.result.folder, 'dir')
    mkdir(param.result.folder)
end

param = setfield(param, 'result', 'logFileName', p.Results.logFileName);
if isempty(param.result.logFileName)
    log_id = getCurrentJob();
    if ~isempty(log_id)
        log_id = ['_', num2str(log_id.ID)];
    else
        log_id = ''; % In case this is not run via batch processing
    end
    [func_stack, ~] = dbstack();
    
    % Log file name
    param.result.logFileName = fullfile(param.result.folder, ...
        sprintf('log_%s_%s.csv', ...
        func_stack(end).file(1:end-2), log_id)); % Log file
else
    param.result.logFileName = fullfile(param.result.folder, ...
        param.result.logFileName); % Log file
end

% Open log file
if isfield(param.result, 'logFileName')
    % Open log file in append mode
    param.result.logFid = fopen(param.result.logFileName, 'a');
    if param.result.logFid == -1
        warning('Cound not open log file')
    end
end

%% Basic link parameters
% Span length
param = setfield(param, 'basicLink', 'spanLength', p.Results.spanLength);

% Dispersion parameter [ps/ns/km], corresponding GVD [s^2/m]
param = setfield(param, 'basicLink', 'D', p.Results.D);
param.basicLink.beta2 = param.basicLink.D*1e6*...
    (-param.constant.wavelength^2/2/pi/param.constant.lightSpeed*1e-6);

% S: Slop of D [s/m^3], correspond to beta3 (TOD) [s^3/m]
param = setfield(param, 'basicLink', 'S', p.Results.S);
param.basicLink.beta3 = (param.basicLink.S-4*pi*param.constant.lightSpeed/...
    param.constant.wavelength^3*param.basicLink.beta2)*...
    (param.constant.wavelength^2/(2*pi*param.constant.lightSpeed))^2*1e-3;

% Attenuation [dB/m]
param = setfield(param, 'basicLink', 'alphaIndB', p.Results.alphaIndB);
% Attenuation [1/m]
param.basicLink.alphaLinear = param.basicLink.alphaIndB*1e3*log(10)/10/1e3;

% Gamma [1/W/m]
param = setfield(param, 'basicLink', 'gamma', p.Results.gamma);

% Number of SSF steps
param = setfield(param, 'basicLink', 'numberSteps', p.Results.numberSteps);

% NF in dB
param = setfield(param, 'basicLink', 'NFdB', p.Results.NFdB);
% Corresponding nsp in linear
param.basicLink.nsp = 10^(param.basicLink.NFdB/10)/2;

% DCF length
param = setfield(param, 'basicLink', 'DCFLength', p.Results.DCFLength);

%% Multiple links
% If linkConfig is an empty cell, setup numberLink links. 
% Else, setup links using the name-value pairs specified in linkConfig
if cellfun(@isempty, p.Results.linkConfig)
    param = setfield(param, 'linkArray', {p.Results.numberLink}, param.basicLink);
else
    linkConfig = p.Results.linkConfig;
    for n=1:length(linkConfig)
        tmp = linkConfig{n};
        nameCell = tmp(1:2:end);
        valueCell = tmp(2:2:end);
        % Use param.basicLink as default values
        param = setfield(param, 'linkArray', {n}, param.basicLink);
        for i=1:length(nameCell)
            param = setfield(param, 'linkArray', {n}, nameCell{i}, valueCell{i});
        end
    end
end

% Number of links
param.numberLink = length(param.linkArray);

%% Basic channel parameters
% Modulation, {'OOK', '16QAM'}
param = setfield(param, 'basicChannel', 'modulation', p.Results.modulation);

% Bandwidth [Hz], equivalent to symbol rate
param = setfield(param, 'basicChannel', 'symbolRate', p.Results.symbolRate);

% Center frequency [Hz]
param = setfield(param, 'basicChannel', 'centerFrequency', ...
    p.Results.centerFrequency);

% Filter parameter
param = setfield(param, 'basicChannel', 'filterFactor', ...
    p.Results.filterFactor);

% Minimum sample per symbol, actual sample per symbol can be larger than
% this number
param = setfield(param, 'basicChannel', 'minSamplePerSymbol', ...
    p.Results.minSamplePerSymbol);

%% Multiple channels
% If channelConfig is an empty cell, setup numberChannel channels with
% default values in basicChannel
% Else, setup channels using the name-value pairs specified in
% channelConfig
if cellfun(@isempty, p.Results.channelConfig)
    param = setfield(param, 'channelArray', {p.Results.numberChannel}, ...
        param.basicChannel);
else
    channelConfig = p.Results.channelConfig;
    for n=1:length(channelConfig)
        tmp = channelConfig{n};
        nameCell = tmp(1:2:end);
        valueCell = tmp(2:2:end);
        % Use param.basicChannel as default values
        param = setfield(param, 'channelArray', {n}, param.basicChannel);
        for i=1:length(nameCell)
            param = setfield(param, 'channelArray', {n}, nameCell{i}, ...
                valueCell{i});
        end
    end
end

% Number of channels
param.numberChannel = length(param.channelArray);

%%
inputs = p.Results;

fprintf(param.result.logFid, 'a log message, %s \n', datestr(now(), 'yyyy-mm-dd HH:MM:SS'));

fclose(param.result.logFid);