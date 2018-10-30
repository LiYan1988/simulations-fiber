function param = getParameters(varargin)
% Set up parameters for simulation.
%
% INPUTS:
%
% OUTPUTS:
%   p: The parameter struct
%
%

%% Parse inputs
% Create an InputParser object
p = inputParser;

checkPositiveScalar = @(x, variableName) assert(isnumeric(x) && ...
    isscalar(x) && (x>0), ...
    sprintf('%s should be positive, scalar, and numeric.', variableName));

checkScalar = @(x, variableName) assert(isnumeric(x) && isscalar(x), ...
    sprintf('%s should be scalar and numeric', variableName));

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
defaultSpanLength = 82; % [km]
checkSpanLength = @(x) checkPositiveScalar(x, 'Span length');
addParameter(p, 'spanLength', defaultSpanLength, checkSpanLength);

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
defaultNumberLinks = 2;
checkNumberLinks = @(x) checkPositiveScalar(x, 'Number of links');
addParameter(p, 'numberLinks', defaultNumberLinks, checkNumberLinks);

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

%% Link parameters
% Number of spans
param = setfield(param, 'link', 'numberSpan', p.Results.numberSpan);

% Span length
param = setfield(param, 'link', 'spanLength', p.Results.spanLength);

% Dispersion parameter [ps/ns/km], corresponding GVD [s^2/m]
param = setfield(param, 'link', 'D', p.Results.D);
param.link.beta2 = param.link.D*1e6*...
    (-param.constant.wavelength^2/2/pi/param.constant.lightSpeed*1e-6);

% S: Slop of D [s/m^3], correspond to beta3 (TOD) [s^3/m]
param = setfield(param, 'link', 'S', p.Results.S);
param.link.beta3 = (param.link.S-4*pi*param.constant.lightSpeed/...
    param.constant.wavelength^3*param.link.beta2)*...
    (param.constant.wavelength^2/(2*pi*param.constant.lightSpeed))^2*1e-3;

% Attenuation [dB/m]
param = setfield(param, 'link', 'alphaIndB', p.Results.alphaIndB);
% Attenuation [1/m]
param.link.alphaLinear = param.link.alphaIndB*1e3*log(10)/10/1e3;

% Gamma [1/W/m]
param = setfield(param, 'link', 'gamma', p.Results.gamma);

% Number of SSF steps
param = setfield(param, 'link', 'numberSteps', p.Results.numberSteps);

% NF in dB
param = setfield(param, 'link', 'NFdB', p.Results.NFdB);
% Corresponding nsp in linear
param.link.nsp = 10^(NF/10)/2;

%% Multiple links
% Number of links
param = setfield(param, 'link', 'numberLinks', p.Results.numberLinks);

% Copy link field to linkArray
param = setfield(param, 'linkArray', {1:param.link.numberLinks}, param.link);

fprintf(param.result.logFid, 'a log message, %s \n', datestr(now(), 'yyyy-mm-dd HH:MM:SS'));

fclose(param.result.logFid);