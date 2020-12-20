% calculate braking period based on heuristics of durations between starts of pulses
% method:
% some pulses are not present, so simple average is not possible.
% - find most common value
% - remove outliers
% - calculate average to get most probable value
% inputs:
% periods - duration between start of pulses, in samples
% ouputs:
% Tb - braking period, in samples

function Tb = get_braking_period(periods)
        m = median(periods);
        % outliers coefficient:
        oc = 0.1;
        periods(periods < m*(1-oc) || periods > m*(1+oc)) = [];
        Tb = mean(periods);
endfunction
