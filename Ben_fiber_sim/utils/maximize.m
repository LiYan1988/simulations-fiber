function maximize(fig)
% Size a window to fill the entire screen.
%
% USAGE:
%    maximize(HANDLE fig)
%
% Modification History
% ??/??/2001  WHF  Created.
% 04/17/2003  WHF  Found 'outerposition' undocumented feature.
%
% Original author: Bill Finger, Creare Inc.
% Free for redistribution as long as credit comments are preserved.

if nargin < 1
    fig = gcf;
end

units = get(fig, 'units');
set(fig, 'units', 'normalized', 'outerposition', [0 0 1 1]);
set(fig, 'units', units);
