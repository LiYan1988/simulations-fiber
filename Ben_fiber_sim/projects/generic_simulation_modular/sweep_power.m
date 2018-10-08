function q = sweep_power()
% MATLAB script to run a generic modularized simulation.
% Script runs a transmitter block, followed by a transmission fiber block,
% which includes a receiver.
%
% Written by Benjamin Foo, 2018-02-19
% Photonics Lab, Chalmers University of Technology

% Signal generation
job_id = 1;
p = sim_transmitter(job_id);

% Creating a parameter sweep
power_array = [15:20];

for sweep_ind = 1:length(power_array)
    p = sim_link(job_id, 'spans', 1, 'power', power_array(sweep_ind));
end

for chan_ind = 1:p.link.N_chan
    q(chan_ind) = sim_postprocess(p.file.path, datestr(p.timestamp, 'yyyy-mm-dd_HH.MM.SS'), chan_ind);
end

end