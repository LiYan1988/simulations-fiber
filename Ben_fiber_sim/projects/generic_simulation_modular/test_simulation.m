function q = test_simulation()
% MATLAB script to run a generic modularized simulation.
% Script runs a transmitter block, followed by a transmission fiber block,
% and then a post-processing block that includes a coherent receiver.
%
% Written by Benjamin Foo, 2018-02-19
% Photonics Lab, Chalmers University of Technology
job_id = 1;
p = sim_transmitter(job_id);

if isfield(p, 'wdm_power') %if there is several WDM powers to simulate, it enters the loop and sets the total wdm power to different values.
    for sweep_ind = 1:length(p.wdm_power)
        p = sim_link(job_id, 'power', p.wdm_power(sweep_ind)); % check how convergence works.
    end
else
    p = sim_link(job_id); 
end

for chan_ind = 6   
    q.pia(chan_ind) = sim_postprocess(p.file.path, datestr(p.timestamp, 'yyyy-mm-dd_HH.MM.SS'), chan_ind, 0);
end
% 
% for chan_ind = 1:p.link.N_chan
%     q.pia(chan_ind) = sim_postprocess(p.file.path, datestr(p.timestamp, 'yyyy-mm-dd_HH.MM.SS'), chan_ind, 0);
% end
% end



% for chan_ind = 1
%     q.pia(chan_ind) = sim_postprocess(p.file.path, datestr(p.timestamp, 'yyyy-mm-dd_HH.MM.SS'), chan_ind, 0);
% end

