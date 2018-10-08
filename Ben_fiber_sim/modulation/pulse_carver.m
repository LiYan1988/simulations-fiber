function u = pulse_carver(p, u)
% Change the incoming NRZ signal to RZ by pulse carving.
%
% INPUTS:
%    p : The parameter struct
%      p.pulse_carver.duty_cycle  : Duty cycle, integer \in {33, 50, 67}
%      p.pulse_carver.unit_energy : Preserve mean power level (boolean)
%    u : The input light signal (1 x M)
% OUTPUTS:
%    u : The modulated light signal (1 x M)
%
% A pulse carver can be driven in three different ways, which will
% determine the duty cycle of the generated pulse. It can be either
% 33%, 50%, or 67%.
%
% See this article for the analytical expressions
%
% @Article{ip_2006_jlt,
%   author = 	 {Ezra Ip and Joseph M. Kahn},
%   title = 	 {Power Spectra of Return-to-Zero Optical Signals},
%   journal = 	 jlt,
%   year = 	 {2006},
%   volume = 	 {24},
%   number = 	 {3},
%   pages = 	 {1610-1617},
%   month = 	 mar,
% }
%
% Pontus Johannisson, 2010-12-16
% This software is distributed under the terms of the GNU General
% Public License version 2

if ~isfield(p.pulse_carver, 'unit_energy');
    % Default is to reduce the signal mean power
    p.pulse_carver.unit_energy = 0;
end;

if p.pulse_carver.unit_energy;
    DC_33 = 0.347878911177953; % See calculation below
    DC_50 = 0.500000000000000;
    DC_66 = 0.652121088822047;
else
    DC_33 = 1;
    DC_50 = 1;
    DC_66 = 1;
end;

t_scaled = p.t/p.dt_symb;
switch p.pulse_carver.duty_cycle;
 case 33;
  envelope = sin(pi/2*(1 + sin(  pi*(t_scaled + 0.5))))/sqrt(DC_33);
 case 50;
  envelope = sin(pi/4*(1 + cos(2*pi*(t_scaled + 0.5))))/sqrt(DC_50);
 case {66, 67};
  envelope = sin(pi/2*(    cos(  pi*(t_scaled + 0.5))))/sqrt(DC_66);
 otherwise;
  error('Undefined duty cycle value');
end;

for row_idx = 1:size(u, 1);
    u(row_idx, :) = u(row_idx, :).*envelope;
end;

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the RZ pulse energy
format long; format compact
t     = linspace(-0.5, 0.5, 1000);
E33   = sin(pi/2*(1 + sin(  pi*(t))));
E50   = sin(pi/4*(1 + cos(2*pi*(t))));
E66   = sin(pi/2*(    cos(  pi*(t))));
DC_33 = trapz(t, E33.^2)
DC_50 = trapz(t, E50.^2)
DC_66 = trapz(t, E66.^2)
