% tries to estimate uncertainty of pulse energy and noise during pulse
% it is based on how correctly the noise around pulse is fitted
% tshift_pulse and tshift_PN are varied based on uncertainty and energy
% of pulse and noise is calculated. output are uncertainties for both
% configurations 1 and 2.

function [uncrE1, uncrE2, uncrEPN1, uncrEPN2] = pulse_uncertainty(Ia, Ib, Vhf, Vdsfhf, fs, ids, ide, tshift_pulse, tshift_PN, tshift_pulse_unc, tshift_PN_unc, plots, groupindex, pulseno, dirpath);

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

end

i = 1; % this is only because I am lazy to delete all '(i)' from equation below, that were copied from energy3.m

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

% calculate relative standard uncertainty:
% the probability distribution function is not simple, so suppose rectangular distribution:
uncrE1   = ( std(E1p(:))./sqrt(3)   )./mean(E1(:));
uncrE2   = ( std(E2p(:))./sqrt(3)   )./mean(E1(:));
uncrEPN1 = ( std(EPN1(:))./sqrt(3) )./mean(E1(:));
uncrEPN2 = ( std(EPN2(:))./sqrt(3) )./mean(E1(:));

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
