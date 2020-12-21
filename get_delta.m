% calculate delta - duty cycle, and its uncertainty
% Inputs:
%       ids - index of sample with pulse start
%       ids - index of sample with pulse end
%       Tb - braking period in samples
% Outputs:
%       delta - duty cycle
%       udelta - uncertainty of duty cycle

function [delta, udelta] = get_delta(ids, ide, Tb)
        delta = (ide - ids)./Tb;
        % estimate of pulse start or end can be by 1 sample
        % off on both start and end, so 2 samples is the
        % uncertainty:
        udelta = 2./Tb;
end
