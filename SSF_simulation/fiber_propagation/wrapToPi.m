%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CS221 Programming Assignment 3
%%   Chris Archibald, Oct. 2009
%%   Stanford University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Wraps an angle so that it is in a consistent format
% This is a function that is used by the simulator internally


function a1 = wrapToPi(a0)
    num2pi = floor(a0/(2*pi) + 0.5);
    a1 = a0 - num2pi*2*pi;
end