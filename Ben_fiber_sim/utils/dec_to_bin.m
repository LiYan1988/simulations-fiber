function [b, n] = dec_to_bin(d, n)
% Convert the input decimal vector to a binary vector 
%
% INPUTS:
%    d : Vector containing decimal values
%    n : Number of bits to use for each decimal value (optional)
% OUTPUTS:
%    b : Vector containing binary values
%    n : Number of bits to use for each decimal value
%
% The size of d is checked, but it is the caller's responsibility to
% make sure that it contains non-negative integers.
%
% Copyright 2011 Pontus Johannisson
% This software is distributed under the terms of the GNU General
% Public License version 2

if (size(d, 1) ~= 1);
    error('Input decimal vector should have size (1 x N)');
end;

if (nargin < 2);
    % Set n so that the largest decimal value can be represented
    max_val = max(d);
    n = floor(log2(max_val)) + 1;
end;

% Set up a look-up table (any better Matlab implementation here???)
b_mat = zeros(n, 2^n);
for k = 1:2^n;
    % Put LSB at highest row index
    b_mat(:, k) = rot90(dec2binvec(k - 1, n), 1);
end;

b = reshape(b_mat(:, d + 1), 1, n*length(d));

return

% Test case
d = [0 1]; n = 1; b = dec_to_bin(d, n) % Should be [0 1]
d = [0 1 2]; n = 2; b = dec_to_bin(d, n) % Should be [0 0 0 1 1 0]
d = [0 1 2 7]; n = 3; b = dec_to_bin(d, n) % Should be [0 0 0 0 0 1 0 1 0 1 1 1]

d = [1 2 3]'; n = 2; b = dec_to_bin(d, n) % Should give error
