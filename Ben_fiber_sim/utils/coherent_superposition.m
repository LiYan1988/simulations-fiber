function [u1_out, u2_out] = coherent_superposition(u1_in, u2_in, test_angle_N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performs a coherent superposition between two optical waveforms.
%
% The function assumes a four-wave mixing-type operation, where one
% waveform is coherently added to the conjugate of the other. The field
% after the coherent superposition is then divided by 2 so that the signal
% power at the input and output should be the same.
%
% For coherent addition, the two waveforms are rotated relative to each
% other by many test angles and the rotation that provides the maximum
% output is selected.
% 
% Inputs:
%   u1_in, u2_in : Optical waveforms to be coherently added
%   test_angle_N : Number of angles to test between 0 and 2*pi rad.
%                  (Default: 3600)
%
% Output:
%   u1_out, u2_out: Optical waveforms after coherent superposition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('test_angle_N', 'var')
    test_angle_N = 3600; % Angular resolution of 0.1 rad
end

if size(u1_in) ~= size(u2_in)
    error('Input wavforms must have the same dimensions')
end

% Testing possibilities for coherent addition
test_angles = linspace(0,2*pi,test_angle_N); % Array for angles to test the coherent superposition
u1_rot_matrix = zeros(1,length(test_angles));
u2_rot_matrix = zeros(1,length(test_angles));

for rot_idx = 1:length(test_angles)
    u1_rot_matrix(rot_idx) = sum(abs(u1_in+conj(u2_in.*exp(1i*test_angles(rot_idx))))); % Sum of field strength? What relevance is the absolute value of the field? You want the result that provides coherent addition (i.e. the optical fields add up in phase)!
    u2_rot_matrix(rot_idx) = sum(abs(u2_in+conj(u1_in.*exp(1i*test_angles(rot_idx)))));
end

% Find maximum of the sum of the field strength - indicates coherent
% ADDITION. Is the index always the same? Should check with more data.
[~,u1_rot_ind] = max(u1_rot_matrix);
[~,u2_rot_ind] = max(u2_rot_matrix);

% Perform the coherent addition (add the fields of the two channels
% in-phase) and reduce field by factor of two so that P_in == P_out
u1_out = (u1_in + conj(u2_in.*exp(1i*test_angles(u1_rot_ind))))./2;
u2_out = (u2_in + conj(u1_in.*exp(1i*test_angles(u2_rot_ind))))./2;

end