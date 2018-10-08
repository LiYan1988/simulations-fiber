function d = bin_to_dec(b, n)
% Convert the input binary matrix to a decimal matrix
%
% INPUTS:
%    b : Matrix containing binary values
%    n : Number of bits to use for each decimal value
% OUTPUTS:
%    d : Matrix containing decimal values
%
% The rows of b are assumed to contain binary data.  These are divided
% into groups of n bits and converted to decimal. This is repeated for
% all rows.
%
% The size of b is checked, but it is the caller's responsibility to
% make sure that it contains binary data. The default value for n is
% to use all bits.
%
% Copyright 2011 Pontus Johannisson
% This software is distributed under the terms of the GNU General
% Public License version 2

if nargin < 2;
    n = size(b, 2);
end;

if (mod(size(b, 2), n))
    error('Input binary vector has incorrect length');
end;

d = zeros(size(b, 1), size(b, 2)/n);
for k = 1:size(b, 1);
    d(k, :) = 2.^(n - 1:-1:0)*reshape(b(k, :), n, size(b, 2)/n);
end;

return

% Test case
b = [0 1]; n = 1; d = bin_to_dec(b, n) % Should be [0 1]
b = [0 0 0 1 1 0]; n = 2; d = bin_to_dec(b, n) % Should be [0 1 2]
b = [0 0 0 0 0 1 0 1 0 1 1 1]; n = 3; d = bin_to_dec(b, n) % Should be [0 1 2 7]
b = [0 0 0 0 0 1 0 1 0 1 1 1; 1 1 1 0 0 0 0 0 1 0 1 0];
n = 3; d = bin_to_dec(b, n) % Should be [0 1 2 7; 7 0 1 2]

b = [0 0 0 1 1 0 0]; n = 2; d = bin_to_dec(b, n) % Should give error
b = [0 0 0 1 1 0]'; n = 2; d = bin_to_dec(b, n) % Should give error
