function raised_cos = raised_cosine(f, dt_symb, beta)
% Generate the frequency response of a raised cosine filter.
% Definition is taken from Wikipedia, but maximum amplitude is one.
%
% INPUTS:
%    f          : Frequency vector [Hz]
%    dt_symbol  : Symbol slot [s]
%    beta       : Roll-off factor (default: 1)
% OUTPUT:
%    raised_cos : Raised cosine frequency response
%
% Pontus Johannisson, 2009-04-29
% This software is distributed under the terms of the GNU General
% Public License version 2

if (nargin < 3)
    beta = 1;
end

if (beta < 0) || (beta > 1)
    error('Roll-off factor outside [0, 1] interval')
end

T = dt_symb;
raised_cos = (abs(f) <= (1 - beta)/2/T);

% Add the roll-off part if beta is finite
if beta > 0
    mask = ((1 - beta)/2/T < abs(f)) .* (abs(f) <= (1 + beta)/2/T);
    val  = 1/2*(1 + cos(pi*T/beta*(abs(f) - (1 - beta)/2/T)));
    raised_cos = raised_cos + val.*mask;
end

%figure(101);clf;stem(f, raised_cos)
