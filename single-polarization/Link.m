classdef Link
    %Fiber link
    %   Detailed explanation goes here
    
    properties (Constant)
        % [m/s], Speed of light
        lightSpeed = 2.99792458e8;
        % [m], reference wavelength
        wavelength = 1550e-9;
        % [J*s] or [W*Hz^-2], Plank's constant
        h = 6.626e-34;
        % [Hz], frequency
        nu = 1.934144890322581e+14;
    end
    
    properties (Access=public)
        % Span length [m]
        spanLength
        % Dispersion parameter [s/m^2]
        D
        % beta2 (GVD) [s^2/m]
        beta2
        % Slop of D [s/m^3]
        S
        % beta3 (TOD) [s^3/m]
        beta3
        % Attenuation [dB/m]
        alphadB
        % Attenuation [1/m]
        alphaLinear
        % Kerr nonlinearity [1/W/m]
        gamma
        % Number of steps
        numberSteps
        % Noise figure in dB [dB]
        NFdB
        % nsp
        nsp
        % DCF length [km]
        DCFLength
    end
    
    methods
        function obj = Link(varargin)
            %Construct an instance of this class
            %   Inputs are name-value pairs, D, S, alphadB, NFdB have
            %   higher priority than beta2, beta3, alphaLinear, nsp
            
            %% Parse input
            p = inputParser;

            addParameter(p, 'spanLength', 82e3, @isnumeric);
            addParameter(p, 'D', 1.7e-5, @isnumeric);
            addParameter(p, 'S', -21.93548387096774, @isnumeric);
            addParameter(p, 'alphadB', 2e-4, @isnumeric);
            addParameter(p, 'gamma', 1.4e-3, @isnumeric);
            addParameter(p, 'NFdB', 5.563, @isnumeric);
            addParameter(p, 'DCFLength', 82e3, @isnumeric);
            addParameter(p, 'numberSteps', 100, @isnumeric);
            
            % Parse the inputs
            parse(p, varargin{:});
            
            %% Set parameters
            % span length
            obj.spanLength = p.Results.spanLength;
            % D and beta2
            obj.D = p.Results.D;
            obj.beta2 = obj.D*1e6*...
                (-obj.wavelength^2/2/pi/obj.lightSpeed*1e-6);
            % S and beta3
            obj.S = p.Results.S;
            obj.beta3 = (obj.S-4*pi*obj.lightSpeed/...
                obj.wavelength^3*obj.beta2)*...
                (obj.wavelength^2/(2*pi*obj.lightSpeed))^2*1e-3;
            % alpha in dB and linear
            obj.alphadB = p.Results.alphadB;
            obj.alphaLinear = obj.alphadB*1e3*log(10)/10/1e3;
            % gamma
            obj.gamma = p.Results.gamma;
            % number of steps in SSF
            obj.numberSteps = p.Results.numberSteps;
            % NF in dB and nsp
            obj.NFdB = p.Results.NFdB;
            obj.nsp = 10^(obj.NFdB/10)/2;
            obj.DCFLength = p.Results.DCFLength;
        end
    end
end

