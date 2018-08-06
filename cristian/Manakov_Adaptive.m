function [Signal, varargout] = Manakov_Adaptive(Signal,P)

% Two axis Adaptive Manakov Equation Solver
% Symmeterised Split step NLSE solver, includes 2nd and 3rd order
% dispersion, PMD, SPM and XGM.
%
% Signal=Manakov_Adaptive(Signal,P)
%
% Inputs:
% SignalIn - input signal structure
% P - Fibre parameters structure contains
%  .Length - fibre length (km)
%  .dz - simulation step size (km)
%  .RefWavelength - reference wavelength (nm)
% Only specify the following parameters if the simulation requires it
%  .Att - fibre attenuation (dB/km)
%  .D - dispersion parameter at reference wavelength (ps/nm/km)
%  .S - dispersion slope at reference wavelength (ps/nm^2/km)
%  .Gamma - nonlinear parameter (/W/km)
%  .PhiMax - Maximum nonlinear phase shift per nonlinear step - rads


% Returns:
% Signal - output signal structure
%
% Author: D Millar - Aug. 2010
% Modified: B Thomsen - Oct. 2010
% Modified: C Czegledi - Sept. 2017
%% IO Verification
P=IOVerify(P);
if P.Length==0
    if isfield(P,'verbose')&&P.verbose>0, disp('Zero length fibre'), end
    P=SItoOpt(-1,P);
    varargout{1} = P;
    return
end

%% Convert SI to Optical Units
P=SItoOpt(1,P);
c=3e5;              % nm/ps


[Np,Nt] = size(Signal.Et); 

FF = MakeTimeFrequencyArray(Signal)/1e12;  % Frequency array [THz]

if (isfield(P, 'Att')&&(P.Att~=0))
    a=P.Att*log(10)/10;         % /km
else
    a=0;
end

if isfield(P,'EDFA')&&~P.EDFA
    a = 0;
    if isfield(P, 'Kt')
        aRamandB = P.Att; % [dB/km]
        aRaman = (aRamandB/4.343);   % [dB/km] --> [Np/km]
        h = 6.62606957e-34; % Kaye and Laby (2015) defines as 6.62606957(29)e-34
        PaseRamanTemp = aRaman*2*h*Signal.Fc*Signal.Fs*P.Kt;     % ASE Power over simulation bandwidth (W);
    else
        aRamandB = 0;
    end    
else
    aRamandB = 0;
end

if isfield(P, 'D')
    B2=-P.D*P.RefWavelength.^2/(2*pi*c);   	                        % ps^2/km
    d2=1i/2*B2*(2*pi*FF).^2;
else
    d2=0;
end

if isfield(P, 'S')
    B3=P.RefWavelength.^2/(2*pi*c).^2*(P.RefWavelength.^2*P.S+2*P.RefWavelength*P.D);      % ps^3/km
    d3=1i*B3/6*(2*pi*FF).^3;
else
    d3=0;
end

clear FF

Po = max(max(abs(Signal.Et).^2));

E=ifft(Signal.Et,[],2);                              % Fourier Transform to frequency domain


if ~isfield(P, 'PhiMax')
    PhiMax = 0.003;
else
    PhiMax = P.PhiMax;
end

if isfield(P, 'Gamma')
    Gamma = P.Gamma;
else
    Gamma = 0;
end

if isfield(P,'verbose')&&P.verbose>0, disp('Using Two Axis Adaptive Split Step NLSE with Manakov Nonlinear estimation'), end


dist=0;
while dist<P.Length
    
    dz=PhiMax/(Gamma*Po);

    if dist+dz>P.Length
        dz=P.Length-dist;
    end
    
    if a>0
        dzEff = (1-exp(-a*dz))/a;
    else
        dzEff = dz;
    end
    
    D = exp(dz/2*(d2+d3));
    E(1,:)=E(1,:).*D;                               % Apply Dispersion operator for the first half step
    E(2,:)=E(2,:).*D;                               % Apply Dispersion operator for the first half step
    E=[fft(E(1,:)); fft(E(2,:))];                            % Fourier Transform to time domain
    
    
%     aEt = abs(Et).^2;
%     NLFactor = 8/9*dzEff*1i*Gamma;
%     N=NLFactor*[aEt(1,:)+aEt(2,:); aEt(2,:)+aEt(1,:)];              % Nonlinear operator SPM & XPM
%     Et=Et.*exp(N-dz*a/2);                                           % Apply Nonlinear and loss operators at center of step

    E=E.*(ones(Np,1)*exp(8/9*dzEff*1i*Gamma*sum(abs(E).^2)-dz*a/2));
    
    if aRamandB>0 % add noise for Raman amplification        
        Pase=dz*PaseRamanTemp;
        % White Gaussian noise zero mean and standard deviation equal to ASE
        % Power split equally across the dimensions
        noiset = sqrt(0.25*Pase)*(randn(Np,Nt)+1i*randn(Np,Nt));
        E = E+noiset;
    end
    
    Po = max(max(abs(E).^2));             % Calculate peak Power
    
    E=[ifft(E(1,:)); ifft(E(2,:))];                    % Fourier Transform to frequency domain    
    E=E.*(ones(Np,1)*D);                               % Apply Dispersion operator for the first half step
    
    dist=dist+dz;
end
Signal.Et=fft(E,[],2);                                                 % Fourier Transform to time domain
% Signal.Et = double(gather(Signal.Et));                                  % GPU --> CPU if required
P=SItoOpt(-1,P);
varargout{1} = P;

end
%% Input/Output Verification
function P=IOVerify(P)
    warningnum = 0;
    if isfield(P,'Att')&&P.Att>1e-3, warning(['Attenuation value is large: ' num2str(P.Att) ' Np/m']), warningnum=warningnum+1; end
    if isfield(P,'D')&&abs(P.D)>1e-3, warning(['Dispersion value is large: ' num2str(P.D) ' s/m^2']), warningnum=warningnum+1; end
    if isfield(P,'Gamma')&&abs(P.Gamma)>1e-1, warning(['Gamma value is large: ' num2str(P.Gamma) ' W/m']), warningnum=warningnum+1; end

    %% TODO
    % PMD, S, dz, Length
    % if P.PMD>1e-1, warning(['PMD value is large: ' num2str(P.PMD) ' s/m^0.5']), warningnum=warningnum+1; end

    %% ENDTODO

    if (P.RefWavelength<1200e-9)||(P.RefWavelength>1700e-9), warning(['wavelength beyond fiber transmission window: ' num2str(P.RefWavelength) ' m']), warningnum=warningnum+1; end

    if (isfield(P,'GPU')&&P.GPU==1&&gpuDeviceCount<1)
        warning('This machine does not support a GPU; Manakov will be computed on the CPU')
        P.GPU=0;
        warningnum=warningnum+1;
    end

    if (warningnum&&isfield(P,'verbose')&&(P.verbose>=1)), fprintf(['<strong>Manakov warnings: ' num2str(warningnum) '\n</strong>']), end
end
