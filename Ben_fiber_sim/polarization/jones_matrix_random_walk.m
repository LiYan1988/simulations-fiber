function U = jones_matrix_random_walk(U, d)
% Perform a random walk to model the evolution of the Jones matrix
%
% INPUTS:
%    U : The Jones matrix from last iteration
%    d : The step parameter
% OUTPUTS:
%    U : The updated Jones matrix
%
% A general Jones matrix can be written
% [ u  v
%  -v* u*],
% with the normalization |u|^2 + |v|^2 = 1.
%
% This implementation aims at being simple and understandable.  Not
% much thought has been given to performance.
%
% Pontus Johannisson, 2009-10-29

% Perform random walk on the sphere in four dimensions
u = U(1, 1);
v = U(1, 2);

% The 4D noise we add has no preferred direction, i.e., it is
% circularly symmeric, i.e., the PDF only depends on the radius.  This
% is important for having a uniform distribution in the long run.
u = u + d*(randn() + 1i*randn());
v = v + d*(randn() + 1i*randn());

% Normalize
n = sqrt(u'*u + v'*v);
u = u/n;
v = v/n;

% Set up updated Jones matrix
U = [ u,  v;
     -v', u'];

% % The following is an alternative model that also should have
% % reasonable properties, see my text about polarizations.  However, we
% % will consider the 4D model above as the main approach.

% % Find a random small rotation
% Omega = d*randn(3, 1);

% % Find the axis/angle
% b_abs = norm(Omega);
% b_hat = Omega/b_abs;

% % Set up the rotation matrix
% T = zeros(2, 2);
% T(1, 1) = cos(b_abs/2)          + i*b_hat(1)*sin(b_abs/2);
% T(1, 2) = b_hat(3)*sin(b_abs/2) + i*b_hat(2)*sin(b_abs/2);
% T(2, 1) = -conj(T(1, 2));
% T(2, 2) =  conj(T(1, 1));

% % Perform the rotation
% U = T*U;

% % For a large number of rotations we can have some error accumulation
% % such that U is not normalized properly any more.  This is not taken
% % into account here.

return

% A simple visualization
close all;
U = eye(2);
d = 0.01;
figure(1); clf;
plot_poincare();
for k = 1:200;
    plot_stokes(jones2stokes(U*[1; 0]));
    U = jones_matrix_random_walk(U, d);
    figure(gcf); drawnow; pause(0.1);
end;
