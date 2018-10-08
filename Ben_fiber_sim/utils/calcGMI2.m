function [GMI,MIperBitPosition]=calcGMI2(X,Y,labeling,N0)
%calcMI calculates the sum of bitwise memoryless mutual informations using circularly symmetric Gaussian noise statistics.
% This quantity is also known as generalized mutual information (GMI). In
% contrast to the symbolwise MI, it is an achievable rate for receivers with binary decoding and no iterations between demapper and decoder.
% The current version of this script works only for uniformly distributed square QAM.
%
% Input:
% X         1 x N       N transmitted complex symbols chosen from 2^m possible values
% Y         1 x N       N received complex symbols
% labeling  chars       optional: 'Gray' (default) or 'Binary'
%
% Output:
% GMI                   1 x 1       Memoryless generalized mutual information assuming circularly symmetric Gaussian noise statistics
% MIperBitPosition      m x 1       Mutual information of each of the m parallel bit channels
%
% Author: Tobias Fehenberger <tobias.fehenberger@tum.de>, Apr. 2015
%
% Downloaded from http://fehenberger.de/#sourcecode
% Modified by Tobias Eriksson <tobias.eriksson@nokia.com>, 20-jun-2016
%

%% number of bits per symbol
m = log2(numel(unique(X)));

%% we want Y as column vector and X as row vector
if size(Y,1) ~= 1
    Y = Y.';
end
if size(X,2) ~= 1
    X = X.';
end

%% we need X also in integer representation
% The input is assumed to be square QAM. For other formats, adapt the demodulator object.
X = X/sqrt(mean(abs(X).^2)); % normalize
hDemod = comm.RectangularQAMDemodulator(2^m, 'BitOutput',false, ...
    'NormalizationMethod', 'Average power', 'SymbolMapping', 'Binary'); % mapping is only important for shaped X when
                                                                        % the symbols must be sorted from top left
                                                                        % running down column-wise to the bottom right.
Xint = step(hDemod,X)';

%% Labeling
if nargin==2
    labeling = 'Gray';
end

%% Calculate sent bits and LLRs
% Y = Y/sqrt(mean(abs(Y).^2));
% N0 = estimateNoiseVar(Xint,Y);

hDemodInput = comm.RectangularQAMDemodulator(2^m, 'BitOutput',true, ...
    'NormalizationMethod', 'Average power', ...
    'SymbolMapping', labeling);

hDemodOutput = comm.RectangularQAMDemodulator(2^m, 'BitOutput',true, ...
    'DecisionMethod', 'Log-likelihood ratio', 'SymbolMapping', labeling,...
    'Variance', N0, 'NormalizationMethod', 'Average power');

c = step(hDemodInput, X)';
LLRs = step(hDemodOutput,Y.')';

%% Compute bitwise MIs and their sum
MIperBitPosition=zeros(m,1);

for kk=1:m
    MIperBitPosition(kk)=1-mean(log2(1+exp((2.*c(kk:m:end)-1).*LLRs(kk:m:end))));
end
GMI=sum(MIperBitPosition);

end

