function u = general_filter(p, u)
% Check the parameter structure and call the specific filtering
% function.
%
% INPUTS:
%    p : The parameter struct
%       p.f                : Frequency vector [Hz]
%       p.filter.type      : Filter function (string)
%       p.filter.bandwidth : Filter bandwidth [Hz]
%       p.filter.bw_spec   : Filter bandwidth spec (string)
%       p.filter.f0        : Center frequency [Hz]
%    u : Input signal [sqrt(W)]
% OUTPUTS:
%    u : Output signal [sqrt(W)]
%
% p.filter.bw_spec can be either 'FWHM' (full width at half maximum)
% or 'HWHM' (half width at half maximum). The latter is more
% convenient for low-pass filters in complex baseband.
%
% Pontus Johannisson, 2010-04-19
% This software is distributed under the terms of the GNU General
% Public License version 2

% Set up default parameters
if ~isfield(p.filter, 'bw_spec')
    p.filter.bw_spec = 'FWHM';
end

if ~isfield(p.filter, 'f0')
    p.filter.f0 = 0;
end

switch p.filter.bw_spec
 case 'FWHM'
  factor_to_FWHM = 1;
 case 'HWHM'
  factor_to_FWHM = 2;
 otherwise
  error('Incorrect p.filter.bw_spec');
end

filt_param.f        = p.f;
filt_param.f0       = p.filter.f0;
filt_param.delta_f  = p.filter.bandwidth*factor_to_FWHM;
switch p.filter.type
 case 'ideal_bp' %Ideal Band Pass Filter
  u = filter_ideal_bp(filt_param, u);
 case 'gaussian' %Gaussian Filter
  u = filter_gaussian(filt_param, u);
 case 'none'
  ;
 otherwise
  error('Incorrect p.filter.type');
end
