%% TODO

%% Things not working
% * *Save generated signals* as |.mat| files for future use.
% 
%       Will not do this. 
%       For Monte Carlo, the data is randomly generated, so no point
%       to save data. 
%       For parameter sweeping, I also need to change parameters and 
%       signals correctly, which also require some effort. 
%
% * Use optimization toolbox in Synchronize to find the optimal sampling
% offset
%
%       Not necessary, optimization works even worse and takes more 
%      function evaluations.

%% Things to try
% * Define functions to plot constellation diagram, eye diagram, SNR/BER
% figures, etc.
%
%  Done
%
% * What to write in log file? Close log file after simulation. 
%
%  Not too much
%
% * Predefine simulation scenarios to make simulation more convenient.
% Or write functions to fast configure simulations.