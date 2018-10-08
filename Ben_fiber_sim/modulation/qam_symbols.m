function qam = qam_symbols()
% Set up normalized N-QAM symbols. (Avoids using the communication
% toolbox.)
%
% INPUTS:
%    none
% OUTPUT:
%    qam : Struct with complex symbols
%
% Pontus Johannisson, 2009-04-29
% This software is distributed under the terms of the GNU General
% Public License version 2

% Mark all places in the code that uses this convention with the
% comment:
% Using conventions from qam_symbols.m
% Should we ever need to change these, it will be valuable to be
% able to find code that depends on the choice made here.

% The follwing definitions are not using the unitary-variance
% assumption. Instead, the symbols have been placed such that the
% constellation is easy to visually inspect.

qam.qam4 = [
    -1 - 1i
    +1 - 1i
    -1 + 1i
    +1 + 1i].';

% With normalized power
qam.qam4_norm = [
    -1 - 1i
    +1 - 1i
    -1 + 1i
    +1 + 1i].'/sqrt(2);

qam.qam16 = [
    -3 - 3i; -1 - 3i; 1 - 3i; 3 - 3i;
    -3 - 1i; -1 - 1i; 1 - 1i; 3 - 1i;
    -3 + 1i; -1 + 1i; 1 + 1i; 3 + 1i;
    -3 + 3i; -1 + 3i; 1 + 3i; 3 + 3i].';

% With normalized power
qam.qam16_norm = [
    -3 - 3i; -1 - 3i; 1 - 3i; 3 - 3i;
    -3 - 1i; -1 - 1i; 1 - 1i; 3 - 1i;
    -3 + 1i; -1 + 1i; 1 + 1i; 3 + 1i;
    -3 + 3i; -1 + 3i; 1 + 3i; 3 + 3i].'/sqrt(10);

% Set up some additional information for the constellations
qam.stat.qam4.Es2 = mean(abs(qam.qam4).^2);
qam.stat.qam4.Es4 = mean(abs(qam.qam4).^4);
qam.stat.qam4.amp = unique(abs(qam.qam4));

qam.stat.qam4_norm.Es2 = mean(abs(qam.qam4_norm).^2);
qam.stat.qam4_norm.Es4 = mean(abs(qam.qam4_norm).^4);
qam.stat.qam4_norm.amp = unique(abs(qam.qam4_norm));

qam.stat.qam16.Es2 = mean(abs(qam.qam16).^2);
qam.stat.qam16.Es4 = mean(abs(qam.qam16).^4);
qam.stat.qam16.amp = unique(abs(qam.qam16));

qam.stat.qam16_norm.Es2 = mean(abs(qam.qam16_norm).^2);
qam.stat.qam16_norm.Es4 = mean(abs(qam.qam16_norm).^4);
qam.stat.qam16_norm.amp = unique(abs(qam.qam16_norm));
