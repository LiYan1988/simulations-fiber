function symnum = data_to_symnum(data)
% Map the data binary matrix to a vector containing the non-binary
% number of the symbols.
%
% INPUTS:
%    data : The binary data (N x M)
% OUTPUTS:
%    symnum : The symbol number (1 x M)
%
% The binary data is 0 or 1. The symbol number is in the interval [1,
% 2^N] since Matlab likes indices to start on 1 and not 0. This means
% the total number of _different_ symbols is 2^N. M is the total
% number of _modulated_ symbols.
%
% P. Johannisson, 2009-05-06
% This software is distributed under the terms of the GNU General
% Public License version 2

N        = size(data, 1);
bin_vect = 2.^((1:N) - 1); % [1 2 4 8 ... 2^N]
symnum   = bin_vect*data + 1;
