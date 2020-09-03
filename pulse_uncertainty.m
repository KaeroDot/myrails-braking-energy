% tries to estimate uncertainty of pulse energy and noise during pulse
% it is based on how correctly the noise around pulse is fitted
% tshift_pulse and tshift_PN are varied based on uncertainty and energy
% of pulse and noise is calculated. output are uncertainties for both
% configurations 1 and 2.

% considered uncertainties:
% - uncertainty for different initial configuration of swticher: Vhf/Ia versus Vdsfhf/Ia
% - uncertainty for detection of pulse start
% - uncertainty for noise level subtraction

function [uncrE1, uncrE2, uncrEPN1, uncrEPN2 report] = pulse_uncertainty(Ia, Ib, Vhf, Vdsfhf, fs, ids, ide, tshift_pulse, tshift_PN, tshift_pulse_unc, tshift_PN_unc, plots, groupindex, pulseno, dirpath);
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
% XXX error of bandwith - from estimation of rectangular peak harmonics cut off
% XXX range? etc.
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

if plots
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
end

i = 1; % this is only because I am lazy to delete all '(i)' from equation below, that were copied from energy3.m

% calculate uncertainty of energy based on uncertainty of noise and pulse limits %<<<1
% generate all possible shifts in time for start/end of pulse
pulse_ids = (tshift_pulse - tshift_pulse_unc):(tshift_pulse + tshift_pulse_unc);
% generate all possible shifts in time for pulse with noise around:
PN_ids = (tshift_PN - tshift_PN_unc):(tshift_PN + tshift_PN_unc);
% loop for all possible starts of pulse:
for k = 1:length(pulse_ids)
        idsS = ids - pulse_ids(k);
        ideS = ide + pulse_ids(k);
        % loop for all possible starts of noise around pulse:
        for j = 1:length(PN_ids)
                idsPN = ids - PN_ids(j);
                idePN = ide + PN_ids(j);
                % copied from energy3 - should be changed into subfunction! XXX

                % copy start here <<<<
                    % x axis - just indexes, sufficient for fit:
                    x = [idsPN(i):idsS(i) ideS(i):idePN(i)];
                    % fit noise around pulse for current a
                    pa = polyfit(x, [Ia(idsPN(i):idsS(i)) Ia(ideS(i):idePN(i))], 1);
                    % tmpIa = Ia(idsS(i):ideS(i)) - polyval(pa, [idsS(i):ideS(i)]);
                    IaN = polyval(pa, [idsS(i):ideS(i)]);
                    % fit noise around pulse for current b
                    pb = polyfit(x, [Ib(idsPN(i):idsS(i)) Ib(ideS(i):idePN(i))], 1);
                    % tmpIb = Ib(idsS(i):ideS(i)) - polyval(pb, [idsS(i):ideS(i)]);
                    IbN = polyval(pb, [idsS(i):ideS(i)]);
                    % energy of pulse:
                    E1(i) = trapz(Vhf   (idsS(i):ideS(i)).*Ia(idsS(i):ideS(i)))./fs + trapz(Vdsfhf(idsS(i):ideS(i)).*Ib(idsS(i):ideS(i)))./fs;
                    E2(i) = trapz(Vdsfhf(idsS(i):ideS(i)).*Ia(idsS(i):ideS(i)))./fs + trapz(Vhf   (idsS(i):ideS(i)).*Ib(idsS(i):ideS(i)))./fs;
                    % energy of noise part during pulse:
                    tmpE1 = trapz(Vhf   (idsS(i):ideS(i)).*IaN)./fs + trapz(Vdsfhf(idsS(i):ideS(i)).*IbN)./fs;
                    tmpE2 = trapz(Vdsfhf(idsS(i):ideS(i)).*IaN)./fs + trapz(Vhf   (idsS(i):ideS(i)).*IbN)./fs;
                    % subtract energy of noise from energy of pulse:
                    E1(i) = E1(i) - tmpE1;
                    E2(i) = E2(i) - tmpE2;
                % copy end here >>>>

                % calculate uncertainty of energy based on current and voltage offsets and gains %<<<1
                % this uncertainty is calculated only for values if idsS,ideS without uncertainties
                % (nonrandomized pulse/noise start/end)
                if (idsS == (ids - tshift_pulse)) && (idsPN == (ids - tshift_PN))
                        % ensure row vector:
                        Iatmp = Ia(idsS(i):ideS(i))(:)';
                        % make monte carlo randomization with offset,gain into matrix (rows are monte carlo)
                        Iam = normrnd(1, urIg, M, 1).*Iatmp + normrnd(0, uIo, M, size(Iatmp, 2));
                        % subtract estimated offset line:
                        Iam = bsxfun(@minus, Iam, IaN);
                        Ibtmp = Ib(idsS(i):ideS(i))(:)';
                        Ibm = normrnd(1, urIg, M, 1).*Ibtmp + normrnd(0, uIo, M, size(Ibtmp, 2));
                        Ibm = bsxfun(@minus, Ibm, IbN);

                        Vhftmp = Vhf(idsS(i):ideS(i))(:)';
                        Vhfm = normrnd(1, urVg, M, 1).*Vhftmp + normrnd(0, uVo, M, size(Vhftmp, 2));
                        Vdsfhftmp = Vdsfhf(idsS(i):ideS(i))(:)';
                        Vdsfhfm = normrnd(1, urVg, M, 1).*Vdsfhftmp + normrnd(0, uVo, M, size(Vdsfhftmp, 2));


                        tmp1 = trapz(Vhfm   .*Iam, 2)./fs;
                        tmp2 = trapz(Vdsfhfm.*Ibm, 2)./fs;
                        uE1 = trapz(Vhfm   .*Iam, 2)./fs + trapz(Vdsfhfm.*Ibm, 2)./fs;
                        uE2 = trapz(Vdsfhfm.*Iam, 2)./fs + trapz(Vhfm   .*Ibm, 2)./fs;
                end
                % put energy into matrix:
                E1p(k,j) = E1(i);
                E2p(k,j) = E2(i);
                EPN1(k,j) = tmpE1;
                EPN2(k,j) = tmpE2;

                if plots
                        plot(haa, idsS(i):ideS(i), IaN)
                        plot(hab, idsS(i):ideS(i), IbN)
                end
        end
end

% calculate relative standard uncertainty: %<<<1
% uncertainty by noise fitting:
% (the probability distribution function is not simple, so suppose rectangular distribution:)
uncrE1   = ( std(E1p(:))./sqrt(3)   )./mean(E1(:));
uncrE2   = ( std(E2p(:))./sqrt(3)   )./mean(E1(:));
uncrEPN1 = ( std(EPN1(:))./sqrt(3) )./mean(E1(:));
uncrEPN2 = ( std(EPN2(:))./sqrt(3) )./mean(E1(:));

% uncertainty caused by offset/gain:
urE1 = (uncrE1.^2 + (std(uE1)./mean(E1))^2 )^0.5;
urE2 = (uncrE2.^2 + (std(uE2)./mean(E2))^2 )^0.5;


uncrE1 = (uncrE1.^2 + urE1.^2)^0.5;
uncrE2 = (uncrE2.^2 + urE2.^2)^0.5;

% display informations: %<<<1
% description is not really clear, improve! XXX
report{end+1} = sprintf('pulse %d. u_r from fitting: %.3g, %.3g, from noise/gain: %.3g, %.3g.', pulseno, uncrE1, uncrE2, urE1, urE2);

if plots
        % finish plotting
        hold off
        set(haa, 'yscale', 'log');
        set(hab, 'yscale', 'log');
        legend(haa, 'Ia', 'pulse start', 'pulse end', 'time shifted pulse start', 'time shifted pulse end', 'start of noise for fit', 'end of noise for fit', 'fit of noise')
        legend(hab, 'Ib', 'pulse start', 'pulse end', 'time shifted pulse start', 'time shifted pulse end', 'start of noise for fit', 'end of noise for fit', 'fit of noise')
        title(haa, sprintf('Gr. %d - Selected current pulse no %d\nNoise in pulse region is approx. %g %%', groupindex, pulseno, mean(EPN1(:))./mean(E1(:)).*100))
        title(hab, sprintf('Gr. %d - Selected current pulse no %d\nNoise in pulse region is approx. %g %%', groupindex, pulseno, mean(EPN2(:))./mean(E2(:)).*100))

        saveplot(sprintf('%05d-selected_current_pulse_Ib_%04d', groupindex, pulseno), dirpath)
        close
        saveplot(sprintf('%05d-selected_current_pulse_Ia_%04d', groupindex, pulseno), dirpath)
        close
end

%% --- Report -------------------- %<<<1
report = strjoin(report, sprintf('\n'));
