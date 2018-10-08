function [GMI]=calcGMI_withNormalization(x,y,labeling)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Caclulates GMI with scaling of the signal power + noise power           %
%                                                                         %
% Input:                                                                  %
% X         1 x N       N transmitted complex symbols chosen from 2^m     %
%                       possible values                                   %
% Y         1 x N       N received complex symbols                        %
% labeling  chars       optional: 'Gray' (default) or 'Binary'            %
%                                                                         %
% Output:                                                                 %
% GMI                   1 x 1       Memoryless generalized mutual         %
%                                   information assuming circularly       %
%                                   symmetric Gaussian noise statistics   %
%                                                                         %
% Written by Tobias Eriksson, tobias.eriksson@nokia.com, 20-jun-2016      %
% Original GMI file obtained from http://fehenberger.de/#sourcecode       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = x/sqrt(mean(abs(x).^2));
h = y*x'/(x*x');
y = y/real(h);
sigma = var(y-x);

% Labeling
if nargin==2
    labeling = 'Gray';
end

GMI = calcGMI2(x,y,labeling,sigma);


end