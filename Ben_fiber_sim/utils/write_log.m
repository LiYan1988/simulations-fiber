function [success] = write_log(p, time_vector, message_str)
% Function to write a time-stamped message to a log file. Should allow most
% functions to maintain compatibility with parameter structures that do not
% have a log file.
%
% Inputs:
%   p          : Parameter structure
%   time_vector: MATLAB time vector for the time stamp - generally now()
%   message_str: String to print to the log file
%
% Outputs:
%   success    : Flag indicating if log was properly written. Returns 1 if 
%                success and 0 if failure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if isfield(p, 'log') && ischar(message_str)
        try
            fprintf(p.log.fid, ['%s, Simulation ID, %s, ', message_str, '\n'], datestr(time_vector, 'yyyy-mm-dd HH:MM:SS'), datestr(p.timestamp, 'yyyy-mm-dd_HH.MM.SS'));
            success = 1;
        catch
            success = 0;
        end
    else
        success = 0;
    end
end