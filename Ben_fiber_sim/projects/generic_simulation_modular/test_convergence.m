function output = test_convergence(steps, power, spans)
% A simple convergence test
% Inputs:
%  steps - an array containing the number of steps per nonlinear length to
%          be tested
%  power - launch power per channel (defaults to value in parameter file)
%  spans - number of spans in the link (defaults to value in parameter
%          file)
%
% Written by Benjamin Foo, 2018-02-19
% Photonics Lab, Chalmers University of Technology

job_id = 1;
p = sim_transmitter(job_id);

if ~exist('power', 'var')
    power = p.link.pwr_per_chan;
end

if ~exist('spans', 'var')
    spans = p.link.N_span;
end

% Sweeping over different step sizes for a given distance and launch power
for sweep_ind = 1:length(steps)
    p = sim_link(job_id, 'spans', spans, 'power', power, 'convergence', steps(sweep_ind));
end

rx_chan = round(p.link.N_chan/2); % Selecting the centre channel
raw_output = sim_postprocess(p.file.path, datestr(p.timestamp, 'yyyy-mm-dd_HH.MM.SS'), rx_chan);

% Reshaping outputs
distance = numel(unique(raw_output.span));
converge_test = numel(steps);
output.nlse_steps = reshape(raw_output.nlse_steps, distance, converge_test);
output.span = reshape(raw_output.span, distance, converge_test);
output.ber = reshape(raw_output.ber, distance, converge_test);
output.evm_db = reshape(raw_output.evm_db, distance, converge_test);
output.gmi = reshape(raw_output.gmi, distance, converge_test);

end