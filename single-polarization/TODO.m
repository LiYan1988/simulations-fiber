% TODO:
% 1. Generate signal and save as .mat files, later simulations can load
% them. No need to regenerate in every simulation. And also save space.
% 2. In parameter struct, use struct of structs to define different types
% of parameters, for example, simulation related, fiber constants, higher
% level control, DSP algorithms, etc.
% 3. Define functions to plot constellation diagram, eye diagram, SNR/BER
% figures, etc.
% 4. Fix bugs in choosing frequency resolution and overall bandwidth.
% 5. Simplify DSP algorithms.
% 6. Higher level configuration of spectrum, channels, etc. 
% 7. Improve efficiency and speed, by reducing unnecessary or redundant
% steps.

% Steps:
% 1. Configure basic/general parameters
% 2. Configure basic simulation environment, like spectrum, time domain
% axis
% 3. Configure specific simulation parameters
% 4. Generate signals
% 5. SSF
% 6. DSP
% 7. Visualize