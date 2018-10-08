function [x, phi] = viterbi_viterbi(p, x)
% This is the Viterbi and Viterbi original algorithm for phase
% estimation
%
% INPUTS:
%    p : The parameter struct
%       p.v_and_v.amp_exp : Amplitude exponent (default: p.v_and_v.rot_symm)
%       p.v_and_v.block_size : Number of symbols per block (default: 64)
%       p.v_and_v.rot_symm : The constellation rotation symmetry (default: 4)
%    x : The samples before phase estimation
% OUTPUTS:
%    x   : The samples after phase estimation
%    phi : The estimated phase
%
% Even if a (2 x N) matrix is input, no joint phase estimation will be
% carried out. The reason is that the polarization demultiplexing
% usually cannot identify the relative phase of the two polarizations.
%
% References:
%
% A. J. Viterbi and A. M. Viterbi, "Nonlinear estimation of
% PSK-modulated carrier phase with application to burst digital
% transmission," IEEE Trans. Inform. Theory, vol. IT-29, no. 4,
% pp. 543-551, July 1983.
%
% B. E. Paden, "A matched nonlinearity for phase estimation of a
% PSK-modulated carrier," IEEE Trans. Inform. Theory, vol. IT-32,
% no. 3, pp. 419-422, May 1986.
%
% Copyright 2011 Pontus Johannisson
% This software is distributed under the terms of the GNU General
% Public License version 2

%disp('ToDo: Handle case when "blocksize" does not divide the number of symbols');

if ~isfield(p, 'v_and_v') || ~isfield(p.v_and_v, 'rot_symm');
    p.v_and_v.rot_symm = 4;
end;

if ~isfield(p, 'v_and_v') || ~isfield(p.v_and_v, 'block_size');
    p.v_and_v.block_size = 64;
end;

if ~isfield(p, 'v_and_v') || ~isfield(p.v_and_v, 'amp_exp');
    p.v_and_v.amp_exp = p.v_and_v.rot_symm;
end;

phi = zeros(size(x));

for k = 1:size(x, 1);
    % Raise to the correct power, reshape, and sum the blocks
    if (p.v_and_v.amp_exp == p.v_and_v.rot_symm)
        tmp = x(k, :).^p.v_and_v.rot_symm;
    else
        % In the paper by Paden (see above), the optimal way to scale the
        % amplitude is investigated
        tmp = abs(x(k, :)).^p.v_and_v.amp_exp.* ...
              exp(i*p.v_and_v.rot_symm*angle(x(k, :)));
    end;

    tmp = sum(reshape(tmp, p.v_and_v.block_size, ...
                      size(x, 2)/p.v_and_v.block_size), 1);

    % Get the phase, unwrap, and divide by the number of levels
    tmp = unwrap(angle(tmp))/p.v_and_v.rot_symm;

    % Duplicate phi values to be able to multiply
    phi(k, :) = kron(tmp, ones(1, p.v_and_v.block_size));

    x(k, :) = x(k, :).*exp(-i*phi(k, :));
end;
