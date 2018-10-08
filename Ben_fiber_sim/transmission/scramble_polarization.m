function [x, A] = scramble_polarization(x)
% Change to a random polarization
%
% INPUTS:
%    x : Input matrix
% OUTPUTS:
%    x : Output matrix
%    A : The Jones matrix used for polarization scrambling
%
% Copyright 2011 Pontus Johannisson
% This software is distributed under the terms of the GNU General
% Public License version 2

% Set up a general scrambling matrix
u = randn() + i*randn();
v = randn() + i*randn();
A = [u,  v; -v', u']/sqrt(u'*u + v'*v);

x = A*x;
