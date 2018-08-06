function [Signal, varargout] = Manakov(Signal, P)
% Two axis Manakov
% Symmeterised Split step Manakov equation solver, includes 2nd and 3rd order
% dispersion, PMD, SPM and XPM.
%
% SignalOut=Manakov(Signal,P)
%
% Inputs:
% Signal - input signal structure
% P - Fibre parameters structure contains
%  .Length - fibre length (m)
%  .dz - simulation step size (m)
%  .RefWavelength - reference wavelength (m)
% Only specify the following parameters if the simulation requires it
%  .Att - fibre attenuation (Np/m)
%  .D - dispersion parameter at reference wavelength (s/m^2)
%  .S - dispersion slope at reference wavelength (s/m^3)
%  .PMD - PMD parameter at reference wavelength (s/m^0.5)
%  .Gamma - nonlinear parameter (/W/m)
%
%  .GPU - [0 1] for GPU computation
%
% Returns:
% SignalOut - output signal structure
%
% Converted from NLSE function by DL (2015) 
%
% Source: C. R. Menyuk, "Nonlinear pulse propagation in birefringent optical
% fibers," IEEE Journal of Quantum Electronics, vol. 23, no. 2, (1987)
%
% see also MANAKOV

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

[Np,Nt] = size(Signal.Et);                  % Total number of points
FF = MakeTimeFrequencyArray(Signal)/1e12;  % Frequency array [THz]


if ~isfield(P, 'Gamma') && ~isfield(P, 'PMD')
    Nz=1;           % single step for linear transmission
    dz=P.Length;
else
    Nz=ceil(P.Length/P.dz);
    dz=P.Length/Nz;
end

if isfield(P, 'Att') && P.Att>0
    a=P.Att*log(10)/10;         % W/km
    % Calculate effective length for nonlinear interaction in presence of fibre loss
    dzEff = (1-exp(-a*dz))/a;
else
    a=0;
    dzEff = dz;
end


if isfield(P,'EDFA')&&~P.EDFA && P.Att>0    
    a = 0;
    if isfield(P, 'Kt')
        aRamandB = P.Att; % [dB/km]
        aRaman = (aRamandB/4.343);   % [dB/km] --> [Np/km]
        h = 6.62606957e-34; % Kaye and Laby (2015) defines as 6.62606957(29)e-34
        PaseRaman = dz*aRaman*2*h*Signal.Fc*Signal.Fs*P.Kt;     % ASE Power over simulation bandwidth (W);
    else
        PaseRaman = 0;
    end      
else
    PaseRaman = 0;
end




if isfield(P, 'D')
    B2=-P.D*P.RefWavelength.^2/(2*pi*c);   	                        % ps^2/km
    d2=1i/2*B2*(2*pi*FF).^2;
else
    d2=zeros(1,size(Signal.Et,2));
end

if isfield(P, 'S')
    B3=P.RefWavelength.^2/(2*pi*c).^2*(P.RefWavelength.^2*P.S+2*P.RefWavelength*P.D);      % ps^3/km
    d3=1i*B3/6*(2*pi*FF).^3;
else
    d3=zeros(1,size(Signal.Et,2));
end

%% Dispersion operator
D = ones(Np,1)*exp(dz/2*(d2+d3));

if Np==1 && isfield(P, 'PMD')
    disp('Single pol with PMD is not implemented! PMD is ignored.')        
elseif isfield(P, 'PMD')
    if isfield(P,'PMDdz') && P.PMDdz>dz
        NzPMD=ceil(P.Length/P.PMDdz);
    else
        NzPMD = Nz;
    end    
    
    PMD_rot_indices = round(linspace(1,Nz+1,NzPMD+1));
    
    DGD_divisions = diff([PMD_rot_indices,Nz+1]);
    DGD_local = 0;

    mean_DGD = P.PMD*sqrt(P.Length);             % mean DGD [ps]
    DGD_dz_mean = mean_DGD/(0.9213*sqrt(NzPMD));    % 0.9213=sqrt(8/(3*pi))
    PMD_sec_ind = 1;
    if isfield(P,'randomDGD') && P.randomDGD
        DGD_dz_std = DGD_dz_mean/5;
    else
        DGD_dz_std = 0;
    end
    DGD_dz = DGD_dz_mean*ones(NzPMD,1)+DGD_dz_std*randn(NzPMD,1); 
    % Rotation coefficients for PMD:
    
    anglesInit = randn(4,NzPMD);
    anglesNormalized = sqrt(sum(anglesInit.^2,1));
    anglesNormalized = anglesInit./repmat(anglesNormalized,size(anglesInit,1),1);
    theta = acos(anglesNormalized(1,:));
    axis = anglesNormalized(2:4,:)./repmat(sin(theta),3,1);
    PMDangles = axis.*repmat(theta,3,1); % these angles scatter the SOP uniformly

    if ~isfield(P,'DGD_dz')
        P.DGD_dz = [];
        P.PMDangles = [];
    end
    P.DGD_dz = [P.DGD_dz; DGD_dz./1e12];
    P.PMDangles = [P.PMDangles, PMDangles];
    
    

    
    
    
    pauliSpins(:,:,1) = [1,  0;...
                         0, -1];
         
    pauliSpins(:,:,2) = [0,  1;...
                         1,  0];
         
    pauliSpins(:,:,3) = [ 0, -1i;...
                         1i,  0];         

    
    rotation = @(angles,pauliSpins) expm(-1i*(angles(1)*pauliSpins(:,:,1)+...
                 angles(2)*pauliSpins(:,:,2)+angles(3)*pauliSpins(:,:,3)));
    
             
                 

        
    
    
    

end

if isfield(P, 'Gamma'),
    Po=max(max(abs(Signal.Et).^2));
    Lnl=1/(P.Gamma*Po);
    if dz > abs(Lnl),
        warning('Step length is longer than the effective nonlinear length');
    end
else
    D=D.^2;
end


if (~isfield(P,'SinglePrecision')||(P.SinglePrecision~=1))
    if (isfield(P,'verbose')&&(P.verbose>=1)), disp('Double precision Manakov (for accuracy)'), end
else
    if (isfield(P,'verbose')&&(P.verbose>=1)), disp('Single precision Manakov (for speed)'), end
        if isfield(P,'PMD')
            PMDangles = single(PMDangles);
            DGD_dz = single(DGD_dz);
        end
    
        Signal.Et=single(Signal.Et);
        D = single(D);
        dz=single(dz);
        dzEff = single(dzEff);
        a = single(a);
end



%% GPU computing section
if isfield(P,'GPU') && (P.GPU==1)
    if isfield(P, 'Gamma') % No point using GPU if not doing FFT
        % Moves DATA on GPU
        if isfield(P,'PMD')
            PMDangles = gpuArray(PMDangles);
            DGD_dz = gpuArray(DGD_dz);
        end
        P.Gamma=gpuArray(P.Gamma);
        Signal.Et=gpuArray(Signal.Et);
        D = gpuArray(D);
        dz=gpuArray(dz);
        dzEff = gpuArray(dzEff);
        a = gpuArray(a);
    end
end


%% Transmission
Ef=ifft(Signal.Et,[],2);                                            % Fourier Transform to frequency domain
if isfield(P, 'PMD') && Np==2
    
    if isfield(P, 'Gamma')
        PMD = [exp(-1i*(DGD_local/(2*2))*2*pi*FF);...
               exp( 1i*(DGD_local/(2*2))*2*pi*FF)];
        if (isfield(P,'verbose')&&(P.verbose>=1)), disp('Using Two Axis Split Step Manakov with SPM, XPM and PMD'), end
        
        NLFactor = 8/9*dzEff*1i*P.Gamma;
        for n=1:Nz
            
            if PMD_rot_indices(PMD_sec_ind)==n
                J = rotation(PMDangles(:,PMD_sec_ind), pauliSpins);             % Make polarisation scrambling matrix
                Ef=J*Ef;                                                        % Apply random polarisation rotation
               
        
                
                DGD_local = DGD_dz(PMD_sec_ind)/DGD_divisions(PMD_sec_ind);
                PMD_sec_ind = PMD_sec_ind + 1;
                PMD = [exp(-1i*(DGD_local/(2*2))*2*pi*FF);...
                       exp( 1i*(DGD_local/(2*2))*2*pi*FF)];              % the second 2 at the denominator comes from the fact that the matrix is applied twice  
            end
            
            Ef = Ef.*D.*PMD;                                            % Apply and Dispersion PMD for the first half            
            
           
%             Et=fft(Ef,[],2);                                            % Fourier Transform to time domain
            Et=[fft(Ef(1,:)); fft(Ef(2,:))];                            % Fourier Transform to time domain
            aEt = abs(Et).^2;
            N=NLFactor*[aEt(1,:)+aEt(2,:); aEt(2,:)+aEt(1,:)];          % Nonlinear operator SPM & XPM
            Et=Et.*exp(N-dz*a/2);                                       % Apply Nonlinear and loss operators at center of step
%             Ef=ifft(Et,[],2);                                           % Fourier Transform to frequency domain

            if PaseRaman>0 % add noise for Raman amplification                    
                % White Gaussian noise zero mean and standard deviation equal to ASE
                % Power split equally across the dimensions
                noiset = sqrt(0.25*PaseRaman)*(randn(Np,Nt)+1i*randn(Np,Nt));
                Et = Et+noiset;
            end
    

            Ef=[ifft(Et(1,:)); ifft(Et(2,:))];                          % Fourier Transform to frequency domain
%             Ef = Ef.*[PMD; 1./PMD];                                     % Apply PMD for the second half
            Ef = Ef.*PMD.*D;                                               % Apply PMD for the second half
            

            

%             Ef=Ef.*D;                                                   % Apply Dispersion for the second half step
        end
    else
        PMD = [exp(-1i*(DGD_local/(2))*2*pi*FF);...
               exp( 1i*(DGD_local/(2))*2*pi*FF)];
        if (isfield(P,'verbose')&&(P.verbose>=1)), disp('Using Two Axis Split Step Manakov with PMD'), end
            for n=1:Nz
                if PMD_rot_indices(PMD_sec_ind)==n
                    J = rotation(PMDangles(:,PMD_sec_ind), pauliSpins);                   % Make polarisation scrambling matrix
                    Ef=J*Ef;                                                    % Apply random polarisation rotation
                    DGD_local = DGD_dz(PMD_sec_ind)/DGD_divisions(PMD_sec_ind);
                    PMD_sec_ind = PMD_sec_ind + 1;
                    PMD = [exp(-1i*(DGD_local/(2))*2*pi*FF);...
                           exp( 1i*(DGD_local/(2))*2*pi*FF)];              % the second 2 at the denominator comes from the fact that the matrix is applied twice  
                end                
                Ef = Ef.*PMD.*D;
                if PaseRaman>0 % add noise for Raman amplification                    
                    % White Gaussian noise zero mean and standard deviation equal to ASE
                    % Power split equally across the dimensions
                    Et=[fft(Ef(1,:)); fft(Ef(2,:))];  
                    noiset = sqrt(0.25*PaseRaman)*(randn(Np,Nt)+1i*randn(Np,Nt));
                    Et = Et+noiset;
                    Ef=[ifft(Et(1,:)); ifft(Et(2,:))];                      
                end
           
            end
    end
elseif isfield(P, 'Gamma') && Np==2
    if (isfield(P,'verbose')&&(P.verbose>=1)), disp('Using Two Axis Split Step Manakov with SPM and XPM'), end
    NLFactor = 8/9*dzEff*1i*P.Gamma;
    for n=1:Nz
        Ef=Ef.*D;                                                       % Apply Dispersion operator for the first half step
%         Et=fft(Ef,[],2);                                                % Fourier Transform to time domain
        Et=[fft(Ef(1,:)); fft(Ef(2,:))];                            % Fourier Transform to time domain
        aEt = abs(Et).^2;
        N=NLFactor*[aEt(1,:)+aEt(2,:); aEt(2,:)+aEt(1,:)];              % Nonlinear operator SPM & XPM
        Et=Et.*exp(N-dz*a/2);                                           % Apply Nonlinear and loss operators at center of step
        
        if PaseRaman>0 % add noise for Raman amplification                    
            % White Gaussian noise zero mean and standard deviation equal to ASE
            % Power split equally across the dimensions
            noiset = sqrt(0.25*PaseRaman)*(randn(Np,Nt)+1i*randn(Np,Nt));
            Et = Et+noiset;
        end
%         Ef=ifft(Et,[],2);                                               % Fourier Transform to frequency domain
        Ef=[ifft(Et(1,:)); ifft(Et(2,:))];                          % Fourier Transform to frequency domain

        Ef=Ef.*D;                                                       % Apply Dispersion operator for the second half step
    end
elseif isfield(P, 'Gamma') && Np==1
    if (isfield(P,'verbose')&&(P.verbose>=1)), disp('Using Single Axis Split Step Manakov'), end
    disp('Manakov: Why are you doing this?')
    NLFactor = dzEff*1i*P.Gamma;
    for n=1:Nz
        Ef=Ef.*D;                                                       % Apply Dispersion operator for the first half step
        Et=fft(Ef,[],2);                                                % Fourier Transform to time domain
%         Et=[fft(Ef(1,:)); fft(Ef(2,:))];                            % Fourier Transform to time domain
        N=9/8*NLFactor*abs(Et).^2;                                      % Nonlinear operator SPM & XPM
        Et=Et.*exp(N-dz*a/2);                                           % Apply Nonlinear and loss operators at center of step
        
        if PaseRaman>0 % add noise for Raman amplification                    
            % White Gaussian noise zero mean and standard deviation equal to ASE
            % Power split equally across the dimensions
            noiset = sqrt(0.25*PaseRaman)*(randn(Np,Nt)+1i*randn(Np,Nt));
            Et = Et+noiset;
        end
        
                
        Ef=ifft(Et,[],2);                                               % Fourier Transform to frequency domain
%         Ef=[ifft(Et(1,:)); ifft(Et(2,:))];                          % Fourier Transform to frequency domain
        Ef=Ef.*D;                                                       % Apply Dispersion operator for the second half step
    end
else
    if (isfield(P,'verbose')&&(P.verbose>=1)), disp('Simulating Dispersion only'), end
    Ef=Ef.*D.*exp(-dz*a/2);                                             % Apply Dispersion and loss operators only
    
    if PaseRaman>0 % add noise for Raman amplification                    
        % White Gaussian noise zero mean and standard deviation equal to ASE
        % Power split equally across the dimensions
        Et=[fft(Ef(1,:)); fft(Ef(2,:))];  
        noiset = sqrt(0.25*PaseRaman)*(randn(Np,Nt)+1i*randn(Np,Nt));
        Et = Et+noiset;
        Ef=[ifft(Et(1,:)); ifft(Et(2,:))];                      
    end
end

Signal.Et=fft(Ef,[],2);                                                 % Fourier Transform to time domain
Signal.Et = double(gather(Signal.Et));                                  % GPU --> CPU if required
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