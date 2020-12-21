% Tries to estimate uncertainty of pulse energy and noise during pulse
% it is based on how correctly the noise around pulse is fitted.
% tshift_pulse and tshift_PN are varied based on uncertainty and energy
% of pulse and noise is calculated. output are uncertainties for both
% configurations 1 and 2.

% considered uncertainties:
% - uncertainty for different initial configuration of swticher: Vhf/Ia versus Vdsfhf/Ia
% - uncertainty for detection of pulse start
% - uncertainty for noise level subtraction

% Script inputs are:
% - Ia: current waveform a
% - Ib: current waveform b
% - Vhf: voltage at half-filter
% - Vdsfhf: voltage after down stream filter
% - fs: data sampling frequency
% - ids: index of pulse start
% - ide: index of pulse end
% - tshift_pulse: time shift (in samples) from pulse start/end to correctly catch whole pulse (safety margin)
% - tshift_PN: time shift (in samples) from pulse start/end to determine noise
% - tshift_pulse_unc: maximum variation of shifts
% - tshift_PN_unc: maximum variation of shifts
% - plots: if set to 1, plotting will happen
% - groupindex: index of a breaking group
% - pulseno: index of a pulse in the breaking group
% - dirpath: directory path to all data files

function [uncrE1, uncrE2, uncrEPN1, uncrEPN2, report] = pulse_uncertainty(Ia, Ib, Vhf, Vdsfhf, fs, delta, udelta, ids, ide, tshift_pulse, tshift_PN, tshift_pulse_unc, tshift_PN_unc, plots, groupindex, pulseno, dirpath);
% uncrE1 - relative uncertainty, configuration 1
% uncrE2 - relative uncertainty, configuration 2
% uncrEPN1 - relative uncertainty of noise, configuration 1
% uncrEPN2 - relative uncertainty of noise, configuration 2
% Ia, Ib - sections of currents
% Vhf, Vdsfhf - sections of voltages
% fs - sampling frequency
% ids - index of pulse start
% ide - index of pulse end
% tshift_pulse - shift of pulse start/end to cover whole pulse
% tshift_PN - shift of pulse start/end to cover noise around for noise level fitting
% tshift_pulse_unc - uncertainty of tshift_pulse
% tshift_PN_unc - uncertainty of tshift_PN
% plots - make plots?
% groupindex - index of group (group of pulses)
% pulseno - pulse number
% dirpath - path for plots

%% CONFIGURATION %%%%%%%%%%%%%%%%%%%%%%%% %<<<1
% current offset uncertainty:
% XXX unknown
uIo = 0.1; % A
% current gain uncertainty:
% from: Monitoring Energy and Power Quality On Board Train, rel. uncertainty 2 %
urIg = 0.02;
% voltage offset uncertainty:
% XXX unknown
uVo = 0.1; % V
% voltage gain uncertainty:
% from: Monitoring Energy and Power Quality On Board Train, rel. uncertainty 0.25 %
urVg = 0.0025;
% Monte Carlo repetitions for gain and offset:
M = 1e4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

report = {};

if plots %<<<1
        % show plot with pulse, a and b current
        idsS = ids - tshift_pulse;
        ideS = ide + tshift_pulse;
        idsPN = ids - tshift_PN;
        idePN = ide + tshift_PN;
        ha = figure;
        % plot will be semilogarithmic, this is to prevent warning messages:
        tmp = Ia;
        tmp(tmp <= 0) = NaN;
        plot([1:length(Ia)], tmp, '-b');
        haa = get(ha, 'children');
        hold on
        hb = figure;
        % plot will be semilogarithmic, this is to prevent warning messages:
        tmp = Ib;
        tmp(tmp <= 0) = NaN;
        plot([1:length(Ib)], tmp, '-b');
        hab = get(hb, 'children');
        hold on
        plot(haa, [ids ids], [1 max(Ia)], '-r');
        plot(haa, [ide ide], [1 max(Ia)], '-r');
        plot(hab, [ids ids], [1 max(Ib)], '-r');
        plot(hab, [ide ide], [1 max(Ib)], '-r');
        plot(haa, [idsS idsS], [1 max(Ia)], '-g');
        plot(haa, [ideS ideS], [1 max(Ia)], '-g');
        plot(hab, [idsS idsS], [1 max(Ib)], '-g');
        plot(hab, [ideS ideS], [1 max(Ib)], '-g');
        plot(haa, [idsPN idsPN], [1 max(Ia)], '-k');
        plot(haa, [idePN idePN], [1 max(Ia)], '-k');
        plot(hab, [idsPN idsPN], [1 max(Ib)], '-k');
        plot(hab, [idePN idePN], [1 max(Ib)], '-k');
        % plot continues ...
end %>>>1

% calculate uncertainty contribution caused by uncertainty of noise and pulse limits %<<<1
% generate all possible shifts in time for start/end of pulse
tshift_pulse_ids = (tshift_pulse - tshift_pulse_unc):(tshift_pulse + tshift_pulse_unc);
% generate all possible shifts in time for pulse with noise around:
tshift_PN_ids = (tshift_PN - tshift_PN_unc):(tshift_PN + tshift_PN_unc);
% loop for all possible starts of pulse:
for k = 1:length(tshift_pulse_ids)
        idsS = ids - tshift_pulse_ids(k);
        ideS = ide + tshift_pulse_ids(k);
        % loop for all possible starts of noise around pulse:
        for j = 1:length(tshift_PN_ids)
                idsPN = ids - tshift_PN_ids(j);
                idePN = ide + tshift_PN_ids(j);
                % RaX, RbX is not needed
                [E1a(k,j), E2a(k,j), EPN1a(k,j), EPN2a(k,j), Ra1, Ra2, Rb1, Rb2, IaN, IbN] = single_pulse_energy(...
                        fs,...
                        Vhf,...
                        Vdsfhf,...
                        Ia,...
                        Ib,...
                        delta,...
                        ideS,...
                        idsS,...
                        idsPN,...
                        idePN,...
                        []);
                if plots
                        % add actual data to plot
                        % ... plot continues
                        plot(haa, idsS:ideS, IaN)
                        plot(hab, idsS:ideS, IbN)
                        % plot continues ...
                end
        end
end

% Calculate uncertainty contribution caused by delta and current and voltage offsets and gains %<<<1
% These uncertainties are calculated only for optimal values of idsS,ideS (nonrandomized pulse/noise
% start/end)
idsS == ids - tshift_pulse;
ideS == ide + tshift_pulse;
idsPN == ids - tshift_PN;
idePN == ide + tshift_PN;

% randomize delta value - uniform distribution:
delta = unifrnd(delta - udelta, delta + udelta, M, 1);
% randomize current gains and offsets
Ia_offset = normrnd(0, uIo, M, size(Ia, 2));
Ia_gain = normrnd(1, urIg, M, 1);
Ib_offset = normrnd(0, uIo, M, size(Ib, 2));
Ib_gain = normrnd(1, urIg, M, 1);
% randomize voltage gains and offsets
Vhf_offset = normrnd(0, uVo, M, size(Vhf, 2));
Vhf_gain = normrnd(1, urVg, M, 1);
Vdsfhf_offset = normrnd(0, uVo, M, size(Vdsfhf, 2));
Vdsfhf_gain = normrnd(1, urVg, M, 1);

% monte carlo loop
for i = 1:M
        % RaX, RbX, IaN, IbN is not needed
        [E1b(i), E2b(i), EPN1b(i), EPN2b(i), Ra1, Ra2, Rb1, Rb2, IaN, IbN] = single_pulse_energy(...
                fs,...
                Vhf.*Vhf_gain(i) + Vhf_offset(i),...
                Vdsfhf.*Vdsfhf_gain(i) + Vdsfhf_offset(i),...
                Ia.*Ia_gain(i) + Ia_offset(i),...
                Ib.*Ib_gain(i) + Ib_offset(i),...
                delta(i),...
                ideS,...
                idsS,...
                idsPN,...
                idePN,...
                []);
end

% calculate total relative standard uncertainty: %<<<1
% uncertainty by noise fitting:
% (the probability distribution function is not simple, so suppose rectangular distribution:)
uncrE1a   = ( std(E1a(:))./sqrt(3)   )./mean(E1a(:));
uncrE2a   = ( std(E2a(:))./sqrt(3)   )./mean(E1a(:));
uncrEPN1a = ( std(EPN1a(:))./sqrt(3) )./mean(EPN1a(:));
uncrEPN2a = ( std(EPN2a(:))./sqrt(3) )./mean(EPN1a(:));

% uncertainty caused by offset/gain:
uncrE1b = std(E1b)./mean(E1b(:));
uncrE2b = std(E2b)./mean(E2b(:));

% total relative uncertainty:
uncrE1 = (uncrE1a.^2 + uncrE1b.^2)^0.5;
uncrE2 = (uncrE2a.^2 + uncrE2b.^2)^0.5;
uncrEPN1 = uncrEPN1a;
uncrEPN2 = uncrEPN2a;

% display informations: %<<<1
% description is not really clear, improve! XXX
report{end+1} = sprintf('pulse %d. u_r from fitting: %.3g, %.3g, from noise/gain: %.3g, %.3g.', pulseno, uncrE1a, uncrE2a, uncrE1b, uncrE2b);
report{end+1} = sprintf('total uncertainty: %.3g, %.3g.', uncrE1, uncrE2);

if plots
        % finish plotting
        hold off
        set(haa, 'yscale', 'log');
        set(hab, 'yscale', 'log');
        legend(haa, 'Ia', 'pulse start', 'pulse end', 'time shifted pulse start', 'time shifted pulse end', 'start of noise for fit', 'end of noise for fit', 'fit of noise')
        legend(hab, 'Ib', 'pulse start', 'pulse end', 'time shifted pulse start', 'time shifted pulse end', 'start of noise for fit', 'end of noise for fit', 'fit of noise')
        title(haa, sprintf('Gr. %d - Selected current pulse no %d\nNoise in pulse region is approx. %g %%', groupindex, pulseno, mean(EPN1a(:))./mean(EPN1a(:)).*100))
        title(hab, sprintf('Gr. %d - Selected current pulse no %d\nNoise in pulse region is approx. %g %%', groupindex, pulseno, mean(EPN2a(:))./mean(EPN1a(:)).*100))

        saveplot(sprintf('%05d-selected_current_pulse_Ib_%04d', groupindex, pulseno), dirpath)
        close
        saveplot(sprintf('%05d-selected_current_pulse_Ia_%04d', groupindex, pulseno), dirpath)
        close
end

%% --- Report -------------------- %<<<1
report = strjoin(report, sprintf('\n'));
