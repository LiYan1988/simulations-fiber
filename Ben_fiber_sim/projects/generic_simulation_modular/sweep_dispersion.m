function q = sweep_dispersion(disp_array, sig_pwr)
% MATLAB script to run a generic modularized simulation.
% Script runs a transmitter block, followed by a transmission fiber block,
% which includes a receiver.
%
% Written by Benjamin Foo, 2018-02-19
% Photonics Lab, Chalmers University of Technology
job_id = 1;
p = sim_transmitter(job_id);

% Creating a parameter sweep
if ~exist('disp_array','var')
    disp_array = [0:5e3:10e3, 11e3:1e3:19e3, 20e3:5e3:70e3]; % Amount of dispersion pre-compensation [m]
end

if ~exist('sig_pwr','var')
    sig_pwr = -10; % Per-channel signal launch power [dBm]
end

for sweep_ind = 1:length(disp_array)
    p = sim_link(job_id, 'spans', 1, 'power', sig_pwr, 'convergence', 1000, 'dispersion', disp_array(sweep_ind));
end

for chan_ind = 1:p.link.N_chan
    q.psa(chan_ind) = sim_postprocess(p.file.path, datestr(p.timestamp, 'yyyy-mm-dd_HH.MM.SS'), chan_ind, 1);
    q.pia(chan_ind) = sim_postprocess(p.file.path, datestr(p.timestamp, 'yyyy-mm-dd_HH.MM.SS'), chan_ind, 0);
end
end