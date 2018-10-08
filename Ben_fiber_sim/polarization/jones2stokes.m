function S = jones2stokes(E)
% Convert input Jones vectors to Stokes vectors
%
% INPUTS:
%    E : Electric field (x and y component) (complex, 2xN)
% OUTPUT:
%    S : The Stokes vectors (real, 3xN)
%
% Pontus Johannisson, 2009-10-14

% This conversion is according to what seems to be the most commonly
% used convention. (For example, the sign of the third component could
% be discussed.)
S = [abs(E(1, :)).^2 - abs(E(2, :)).^2;
     2*real(E(1, :).*conj(E(2, :)));
    -2*imag(E(1, :).*conj(E(2, :)))];

return
% Test cases, using "Polarization Optics in Telecommunications" by Damask, p. 35
jones2stokes([1;  0])         % [ 1;  0;  0], linear, horizontal
jones2stokes([0;  1])         % [-1;  0;  0], linear, vertical
jones2stokes([1;  1]/sqrt(2)) % [ 0;  1;  0], linear, +45
jones2stokes([1; -1]/sqrt(2)) % [ 0; -1;  0], linear, -45
jones2stokes([1;  i]/sqrt(2)) % [ 0;  0;  1], circular, right
jones2stokes([1; -i]/sqrt(2)) % [ 0;  0; -1], circular, left
% Given the convention above, the mapping between the Jones and Stokes
% vectors is given. However, the designation "right" and "left" for
% the circular polarization states is a further non-universal
% convention. For example, "Optics" by Hecht uses the opposite choice.

% A further test case
jones2stokes([sqrt(2)  0       1  1  1  1;
              0        sqrt(2) 1 -1  i -i]/sqrt(2))
