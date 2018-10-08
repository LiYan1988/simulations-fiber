function e = electric_signal_rrcf(p,k)
% Generate one electric signal.
%
% INPUTS:
%    p : The parameter struct
%       p.electric.data   : Binary data (R x S)
%       p.electric.levels : Electric modulation levels (1 x 2^R)
% OUTPUTS:
%    e : The generated electric signal (1 x M)
%
% The electric signal source has R binary input streams and 1 output
% signal. All of the R binary streams contribute one bit per symbol.
% Modulation is done between the levels in p.electric.levels as
% selected by the data streams.
%
% In order to obtain an electric signal with the correct spectrum, the
% electrical levels are first used to make an ideal wide-band signal,
% and then a lowpass filter is applied.
%
% Pontus Johannisson, 2010-04-19
% This software is distributed under the terms of the GNU General
% Public License version 2

level = data_to_symnum(p.electric.data);

% Ideal signal, one sample per symbol slot
e = p.electric.levels(level);

% Duplicate the ideal signals to the sampling rate
%e = kron(e, ones(1, p.samp_per_symbol(k).wdm));




% Shift the electric signal half a sample duration to get it centered
% exactly in the middle of the symbol slot
e = 0.5*(e + circshift(e, [0 1]));

%Perform the specified filtering
% p.filter = p.electric.filter;
% e = general_filter(p, e);
