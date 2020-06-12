%% === Calculate energy of a braking group ==================== %<<<1
% A rise and fall of pulses of first current waveform are found. Time shift of few samples before and
% after is added and at these points the current waveforms are split. More time shift is added to
% estimate current offset around the pulse.
% The energy is calculated as integral in a whole part between these two split points for pulses
% minus the current offset.
% The energy where pulse is not present is accounted as energy of "noise".
% The two voltages are switched for each pulse for two currents. For the first pulse the first
% voltage is applied to first current -> configuration 1. For the first pulse the second voltage is
% applied to first current -> configuration 2. 
function [E EPN EN urE] = energy3(groupindex, fs, triglvl, dirpath, plots);

%% CONFIGURATION %%%%%%%%%%%%%%%%%%%%%%%% %<<<1
% Time shift by which will be shifted back in time the point of split for start of pulse and shifted forward the
% point of split of end of pulse. Pulse is found by middle of peak - but in real it starts few
% samples back.
tshift_pulse = 20; % (samples) 
tshift_pulse_unc = 5; % (samples) 
% Maximum time shift that will be used to find noise level around the peak. If distance between
% peaks is too small, the time shift will be set to smaller number.
tshift_PN = 50; % (samples) 
tshift_PN_unc = 5; % (samples) 
% how many pulses will be evaluated for uncertainty (and plotted if plots=1):
noofpulses_for_unc = 100;
% randomize which pulses will be evaluated for uncertainty (and plotted if plots=1)?:
% (set 0 for debugging)
randomize_pulses_for_unc = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% --- Load data -------------------- %<<<1
varnms = {'Ia', 'Ib', 'Vdsf', 'Vhf'};
% generate all filenames
for j = 1:length(varnms)
        fng{j} = fullfile(dirpath, [varnms{j} sprintf('-%05d', groupindex) '.mat']);
        load(fng{j});
end

%% --- Other required quantities -------------------- %<<<1
% time axis:
t = [1:length(Ia)]./fs - 1/fs;
% half-filter voltage:
% complementary half-filter voltage:
Vdsfhf = Vdsf - Vhf;

%% --- Display basic informations about signal -------------------- %<<<1
disp(['Length of signal: ' num2str(length(Ia)) ' samples, duration: ' num2str(length(Ia)./fs) ' s.'])
disp(sprintf('Ia: max: %g A, std: %g A.', max(Ia), std(Ia)))
disp(sprintf('Ib: max: %g A, std: %g A.', max(Ib), std(Ib)))
disp(sprintf('Vhf: mean: %g V, std: %g V.', mean(Vhf), std(Vhf)))
disp(sprintf('Vdsf: mean: %g V, std: %g V.', mean(Vdsf), std(Vdsf)))

%% --- Finding indexes of pulse starts and pulse ends and pulse lengths -------------------- %<<<1
% all below trigger is zero, all above is 1:
braking = Ia > triglvl; % contains only 0, 1
% calculate derivation to find out where signal rise up and fall down:
slopes = diff(braking); % contains only -1, 0, 1

% get indexes of all slopes:
allslopes = slopes;
allslopes(allslopes == -1) = 1;
idslopes = find(allslopes); 

% get indexes of positive slopes (rise up):
posslopes = slopes;
posslopes(posslopes == -1) = 0;
ids = find(posslopes);

% get indexes of negative slopes (rise down):
negslopes = slopes;
negslopes(negslopes == 1) = 0;
ide = find(negslopes);

% check the same number of rise ups and rise downs:
if length(ids) ~= length(ide)
        error('The number of pulse starts is different to number of pulse ends. The splitting into groups was not done correctly.')
end
pulses = 0;
if ~isempty(ids)
        % boolean if pulses are present
        pulses = 1;
end

% get lengths between all slopes (in samples): 
lengths = diff(idslopes);

% get lengths between positive and negative slopes - when voltage is high (in samples): 
tmp = ismember(idslopes, ids);
lengthson = lengths(tmp(1:end-1));

% get lengths between negative and positive slopes - when voltage is low (in samples): 
tmp = ismember(idslopes, ide);
lengthsoff = lengths(tmp(1:end-1));

% get lengths between positive slopes - periods of pulses (in samples): 
periods = diff(ids);

%% --- Display informations about pulses -------------------- %<<<1
if pulses
        disp('Current Ia:')
        disp(['off lengths in samples: ' ...
                  'count=' num2str(length(lengthsoff)) ...
                ', mean='  num2str(mean(lengthsoff)) ...
                ', median=' num2str(median(lengthsoff)) ...
                ', std=' num2str(std(lengthsoff)) ...
                ', min=' num2str(min(lengthsoff)) ...
                ', max=' num2str(max(lengthsoff))])
        disp(['on lengths in samples: ' ...
                  'count=' num2str(length(lengthson)) ...
                ', mean='  num2str(mean(lengthson)) ...
                ', median=' num2str(median(lengthson)) ...
                ', std=' num2str(std(lengthson)) ...
                ', min=' num2str(min(lengthson)) ...
                ', max=' num2str(max(lengthson))])
        disp(['periods in samples: ' ...
                  'count=' num2str(length(periods)) ...
                ', mean='  num2str(mean(periods)) ...
                ', median=' num2str(median(periods)) ...
                ', std=' num2str(std(periods)) ...
                ', min=' num2str(min(periods)) ...
                ', max=' num2str(max(periods))])
        disp(['off lengths in miliseconds: ' ...
                  'count=' num2str(length(lengthsoff)) ...
                ', mean='  num2str(mean(lengthsoff)./fs*1000) ...
                ', median=' num2str(median(lengthsoff)./fs*1000) ...
                ', std=' num2str(std(lengthsoff./fs*1000)) ...
                ', min=' num2str(min(lengthsoff./fs*1000)) ...
                ', max=' num2str(max(lengthsoff./fs*1000))])
        disp(['on lengths in miliseconds: ' ...
                  'count=' num2str(length(lengthson)) ...
                ', mean='  num2str(mean(lengthson)./fs*1000) ...
                ', median=' num2str(median(lengthson)./fs*1000) ...
                ', std=' num2str(std(lengthson)./fs*1000) ...
                ', min=' num2str(min(lengthson)./fs*1000) ...
                ', max=' num2str(max(lengthson)./fs*1000)])
        disp(['periods in miliseconds: ' ...
                  'count=' num2str(length(periods)) ...
                ', mean='  num2str(mean(periods)./fs*1000) ...
                ', median=' num2str(median(periods)./fs*1000) ...
                ', std=' num2str(std(periods)./fs*1000) ...
                ', min=' num2str(min(periods)./fs*1000) ...
                ', max=' num2str(max(periods)./fs*1000)])
else
        disp('No detected pulses. Only noise data');
end

%% --- Get split indexes and times -------------------- %<<<1
% description of indexes of pulses and noise sections:
%     <--> lengthson - durations of on state
%         <--------> lengthsoff - durations of off (noise) state
%               <>        <> tshift_PN
%                 <>    <> tshift_pulse
%      __            __
%     |  |          |  |
%   __|  |__________|  |____
%               ^ ^ ^  ^ ^ ^
%               | | |  | | |
%               | | |  | |  idePN - end of section for noise level estimation
%               | | |  |  ideS - end of pulse shifted by tshift_pulse. also start of section for noise level estimation
%               | | |   ide - end of pulse, really midle of pulse falling slope
%               | |  ids - start of pulse, really middle of pulse rising slope
%               |  idsS - start of pulse shifted by tshift_pulse. also end of section for noise level estimation
%                idsPN - start of section for noise level estimation

if pulses
        % pulse/PN shift indexes %<<<2
        tshift_pulse = min([tshift_pulse fix((lengthsoff - 1)./2)]);
        tshift_PN = min([tshift_PN fix((lengthsoff - 1)./2)]);
        disp(['Minimal off length is ' num2str(min(lengthsoff)) ' samples. All pulse starts will be shifted back by ' num2str(tshift_pulse) ' samples. All pulse ends will be shifted forward by ' num2str(tshift_pulse) ' samples.'])
        disp(['The noise around pulse will be estimated from ' num2str(tshift_PN) ' samples before and after pulse start.'])

        % pulse start indexes %<<<2
        idsS = ids - tshift_pulse;
        % fix first negative index:
        if idsS(1) < 1
                idsS(1) = 1;
        end
        ideS = ide + tshift_pulse;
        % fix last too big index:
        if ideS(end) > length(slopes)
                ideS(end) = length(slopes);
        end

        % pulse noise indexes %<<<2
        idsPN = idsS - tshift_PN;
        % fix first negative index:
        if idsPN(1) < 1
                idsPN(1) = 1;
        end
        idePN = ideS + tshift_PN;

        % get time of splittings %<<<2
        t_pulsestart = t(ids);
        t_pulseend = t(ide);
        t_pulsestartO = t(idsS);
        t_pulseendO = t(ideS);

        % fit time starts by line
        [P, S] = polyfit ([1:length(t_pulsestart)], t_pulsestart, 1);
else
        % no pulses, no splitting %<<<2
        ids = length(Ia);
        ide = ids;
        idsS = ids;
        ideS = ids;
        idsPN = ids;
        idePN = ids;
end

%% --- make calculation in a loop for every pulse -------------------- %<<<1
% following for loop could have been done by a series of arrayfun and cellfun, however it takes 20 % longer time
% than one single for loop

% preparation for plotting, for selected pulses fit data will be saved:
if pulses
        if randomize_pulses_for_unc
                % randomly generate which pulses will be selected for, non repeating indexes:
                pulses_for_unc = ids(randperm(length(ids)))(1:noofpulses_for_unc);
        else
                pulses_for_unc = fix(linspace(2, length(ids)-1, noofpulses_for_unc));
        end
end
% prepare variable to speed up for loop:
E1 = NaN.*zeros(1, size(ids,2));
uncrE1 = E1;
E2 = E1;
uncrE2 = E1;
EN1 = E1;
EN2 = E1;
uncrEPN1 = E1;
uncrEPN2 = E1;
% noise of the first section:
if idsS(1) > 1
        % there is some noise before first pulse:
        EN1(1) = trapz(Vhf   (1:idsS(1)).*Ia(1:idsS(1)))./fs + trapz(Vdsfhf(1:idsS(1)).*Ib(1:idsS(1)))./fs;
        EN2(1) = trapz(Vdsfhf(1:idsS(1)).*Ia(1:idsS(1)))./fs + trapz(Vhf   (1:idsS(1)).*Ib(1:idsS(1)))./fs;
end
if pulses
        % main for loop
        % add start of next fictive pulse as last index to calculate in also last section of noise
        % if idsS(end) is already equal to length(Ia), nothing happen because there will be only 0 value.
        idsS(end+1) = length(Ia);
        for i = 1:(length(ids))
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
                % resistance of shunts during pulse:
                id = find(max(Ia(idsS(i):ideS(i))) == Ia(idsS(i):ideS(i)));
                id = id(1) + idsS(i) - 1;
                Ra1(i) = Vhf   (id)./Ia(id);
                Ra2(i) = Vdsfhf(id)./Ia(id);
                Rb1(i) = Vhf   (id)./Ib(id);
                Rb2(i) = Vdsfhf(id)./Ib(id);
                % energy of pulse:
                E1(i) = trapz(Vhf   (idsS(i):ideS(i)).*Ia(idsS(i):ideS(i)))./fs + trapz(Vdsfhf(idsS(i):ideS(i)).*Ib(idsS(i):ideS(i)))./fs;
                E2(i) = trapz(Vdsfhf(idsS(i):ideS(i)).*Ia(idsS(i):ideS(i)))./fs + trapz(Vhf   (idsS(i):ideS(i)).*Ib(idsS(i):ideS(i)))./fs;
                % energy of noise part during pulse:
                EPN1(i) = trapz(Vhf   (idsS(i):ideS(i)).*IaN)./fs + trapz(Vdsfhf(idsS(i):ideS(i)).*IbN)./fs;
                EPN2(i) = trapz(Vdsfhf(idsS(i):ideS(i)).*IaN)./fs + trapz(Vhf   (idsS(i):ideS(i)).*IbN)./fs;
                % subtract energy of noise from energy of pulse:
                E1(i) = E1(i) - EPN1(i);
                E2(i) = E2(i) - EPN2(i);
                % energy of noise after pulse:
                EN1(i+1) = trapz(Vhf   (ideS(i):idsS(i+1)).*Ia(ideS(i):idsS(i+1)))./fs + trapz(Vdsfhf(ideS(i):idsS(i+1)).*Ib(ideS(i):idsS(i+1)))./fs;
                EN2(i+1) = trapz(Vdsfhf(ideS(i):idsS(i+1)).*Ia(ideS(i):idsS(i+1)))./fs + trapz(Vhf   (ideS(i):idsS(i+1)).*Ib(ideS(i):idsS(i+1)))./fs;
                % add energy of noise during pulse into energy of noise after pulse:
                EN1(i+1) = EN1(i+1) + EPN1(i);
                EN2(i+1) = EN2(i+1) + EPN2(i);
                % uncertainty of estimation
                % uncertainty of estimation does the plotting of individual pulses
                % cut part of currents/voltages so small data are transfered into pulse_uncertainty function:
                idx1 = idsPN(i) - 2*tshift_PN_unc;
                idx2 = idePN(i) + 2*tshift_PN_unc;
                if any(ids(i) == pulses_for_unc)
                        [uncrE1(i), uncrE2(i), uncrEPN1(i), uncrEPN2(i)] = pulse_uncertainty(Ia(idx1:idx2), Ib(idx1:idx2), Vhf(idx1:idx2), Vdsfhf(idx1:idx2), fs, ids(i) - idx1 + 1, ide(i) - idx1 + 1, tshift_pulse, tshift_PN, tshift_pulse_unc, tshift_PN_unc, plots, groupindex, i, dirpath);
                end
        end
        % remove start of fictive pulse that was added before the for loop
        idsS(end) = [];
end % if pulses

%% --- total energy -------------------- %<<<1
% calculate total energy for both configurations (a*hf, b*dsfhf):
if pulses
        % prepare indexing vectors for two powers 
        % i.e. make series 1,0,1,0,... (based on square.m from statistics package):
        sw1 = [1:length(E1)];
        sw1 = 0.5*sw1;
        sw1 = (sw1 - floor(sw1)>=0.5);
        E(1,1) = sum(     sw1 .*E1  + not(sw1).*E2 );
        E(2,1) = sum( not(sw1).*E1 +      sw1 .*E2 );
        EPN(1,1) = sum(     sw1 .*EPN1  + not(sw1).*EPN2 );
        EPN(2,1) = sum( not(sw1).*EPN1 +      sw1 .*EPN2 );
        % resistors:
        Ra(1,:) =     sw1 .*Ra1  + not(sw1).*Ra2;
        Ra(2,:) = not(sw1).*Ra1  +     sw1 .*Ra2;
        Rb(1,:) =     sw1 .*Rb1  + not(sw1).*Rb2;
        Rb(2,:) = not(sw1).*Rb1  +     sw1 .*Rb2;
else
        E(1,1) = 0;
        E(2,1) = 0;
        EPN(1,1) = 0;
        EPN(2,1) = 0;
end
% prepare indexing vectors for two powers 
% i.e. make series 1,0,1,0,... (based on square.m from statistics package):
sw1N = [1:length(EN1)];
sw1N = 0.5*sw1N;
sw1N = (sw1N - floor(sw1N)>=0.5);
EN(1,1) = sum(     sw1N .*EN1  + not(sw1N).*EN2 );
EN(2,1) = sum( not(sw1N).*EN1 +      sw1N .*EN2 );

% total uncertainties: %<<<1
% (lengthy code because isnan is not in basic matlab)
if pulses
        % absolute uncertainties of pulses:
        % simple sum (not sqrt of sum of squares) because correlation is considered
        % as 1, therefore sum of input quantities propagates into simple sum of uncertainties
        % (see GUM Guide 1, page 21 top)
        uE(1,1) = nansum(     sw1 .*uncrE1.*E1 + not(sw1).*uncrE2.*E2);
        uE(2,1) = nansum( not(sw1).*uncrE1.*E1 +     sw1 .*uncrE2.*E2);
        % multiply by ratio of lengths to get approximate uncertainty for all
        % pulses (whole energy), not only the ones randomly selected for uncertainty
        % calculation
        uE(1,1) = uE(1).*length(ids)./length(pulses_for_unc);
        uE(2,1) = uE(2).*length(ids)./length(pulses_for_unc);
        urE = uE./E;
else
        urE = [0; 0];
end

%% --- Display info about energies -------------------- %<<<1

disp(['Energy of total noise (overestimated) conf. 1: ' num2str(EN(1)) ' J, conf. 2: ' num2str(EN(2)) ' J.']);
disp(['Error between total noise energies of two configurations (EN2-EN1)/EN2 (%): ' num2str((EN(2)-EN(1))/EN(1).*100)]);
disp(['Error between total noise energies of two configurations (EN1-EN2)/EN1 (%): ' num2str((EN(1)-EN(2))/EN(2).*100)]);
if pulses
        disp(['Energy of noise during pulse conf. 1: ' num2str(EPN(1)) ' J, conf. 2: ' num2str(EPN(2)) ' J.']);
        disp(['Ratio of noise during pulse to pulse energy, conf. 1 is ' num2str(sum(EPN(1))./sum(E(1)).*100) ' %.']);
        disp(['Ratio of noise during pulse to pulse energy, conf. 2 is ' num2str(sum(EPN(2))./sum(E(2)).*100) ' %.']);

        disp(['Total noise conf. 1 is ' num2str(sum(EN(1))./sum(E(1)).*100) ' % of total pulse energy conf. 1.']);
        disp(['Total noise conf. 2 is ' num2str(sum(EN(2))./sum(E(2)).*100) ' % of total pulse energy conf. 2.']);

        disp(['Error between total pulse energies of two configurations (E2-E1)/E2 (%): ' num2str((E(2)-E(1))/E(1).*100)]);
        disp(['Error between total pulse energies of two configurations (E1-E2)/E1 (%): ' num2str((E(1)-E(2))/E(2).*100)]);

        disp(['Relative uncertainty of single pulse energy E1: mean: ' num2str(nanmean(uncrE1).*100) ', min: ' num2str(min(uncrE1).*100) ', max: ' num2str(max(uncrE1).*100) ' %']);
        disp(['Relative uncertainty of single pulse energy E2  mean: ' num2str(nanmean(uncrE2).*100) ', min: ' num2str(min(uncrE2).*100) ', max: ' num2str(max(uncrE2).*100) ' %']);

        disp(['>> Total pulse energy conf. 1: ' num2str(E(1)) ' J, conf. 2: ' num2str(E(2)) ' J.']);
        disp(['>> Uncertainty of total pulse energy: conf. 1: ' num2str(uE(1)) ' J, conf. 2: ' num2str(min(uE(2))) 'J, (' num2str(uE(1)./E(1).*100) ' %, ' num2str(uE(2)./E(2).*100) ' %).']);

        Esim(1,1) = sum(lengthson)./fs.*max(Ia).*mean(Vhf) + sum(lengthson)./fs.*max(Ib).*mean(Vdsfhf);
        Esim(2,1) = sum(lengthson)./fs.*max(Ia).*mean(Vdsfhf) + sum(lengthson)./fs.*max(Ib).*mean(Vhf);
        disp(['Total pulse energy based on pulse length and average I,V, for Ia, Vhf: ' num2str(Esim(1)) ' J, for Ib, Vdsfhf: ' num2str(Esim(2)) ' J,'])
        disp(['first value is ' num2str(Esim(1)./E(1).*100) ' % of total pulse energy conf. 1.'])
end

%% --- 2DO -------------------- XXX %<<<1
% calculate average offset level?

%% --- Plotting -------------------- %<<<1
if plots
        % plot current with lines showing splitting  %<<<2
        % this figure is challenging the hardware, use only if needed!
        figure
        hold on
        if pulses
                yl = [min(Ia) max(Ia)];
                yl(1) = 1;
                disp('this will be slow...')
                y = repmat(yl', 1, length(t(ids)));
                x = repmat(t(ids), 2, 1);
                plot(x,y,'-r');
                x = repmat(t(ide), 2, 1);
                disp('still working...')
                plot(x,y,'-r');
                x = repmat(t(idsS), 2, 1);
                disp('still working....')
                plot(x,y,'-g');
                x = repmat(t(ideS), 2, 1);
                disp('still working.....')
                plot(x,y,'-g');
                x = repmat(t(idsPN), 2, 1);
                disp('still working......')
                plot(x,y,'+k');
                x = repmat(t(idePN), 2, 1);
                disp('still working.......')
                plot(x,y,'ok');
                title(sprintf('Gr. %d - Current waveform Ia\nmean: %g A, std: %g A\nred: pulse, green: shifted pulse, black +o: noise for pulse', groupindex, mean(Ia), std(Ia)))
        else
                title(sprintf('Gr. %d - Current waveform Ia\nmean: %g A, std: %g A', groupindex, mean(Ia), std(Ia)))
        end
        plot(t, Ia, '-b')
        hold off
        xlabel('time (s)')
        ylabel('I (A)')
        saveplot(sprintf('%05d-current_Ia', groupindex), dirpath)
        close

        if pulses
                % errors of pulse starts from ideal line %<<<2
                % reveals gaps between pulses
                figure
                hold on
                plot(polyval(P, [1:length(t_pulsestart)]) - t_pulsestart, '-b')
                plot(xlim, 0 + 0.5*[min(lengthsoff) min(lengthsoff)]./fs, 'r--');
                plot(xlim, 0 - 0.5*[min(lengthsoff) min(lengthsoff)]./fs, 'r--');
                legend('error', '0 - 1/2 of minimal "off" length', '0 - 1/2 of minimal "off" length', 'location', 'southwest')
                title(sprintf('Gr. %d - Errors of pulse starts time from ideal line fit\nminimal off length: %g samples', groupindex, min(lengthsoff)))
                xlabel('pulse start (s)')
                ylabel('error from line fit (s)')
                saveplot(sprintf('%05d-pulse_start_times-errors', groupindex), dirpath)
                close

                % % energy in splits - pulses
                % figure
                % x = t_pulsestart;
                % y = Ehf_a.*sw1;
                % x(y == 0) = [];
                % y(y == 0) = [];
                % x2 = t_pulsestart;
                % y2 = Edsfhf_a.*not(sw1);
                % x2(y2 == 0) = [];
                % y2(y2 == 0) = [];
                % plot(x, y, 'xr', x2, y2, 'xb')
                % legend('Ehf\_a', 'Edsfhf\_a')
                % title([num2str(groupindex) ' - Energy of pulses, current a, configuration 1'])
                % xlabel('time of pulse start (s)')
                % ylabel('Energy (J)')
                % saveplot(sprintf('%05d-energy_a_config_1', groupindex), dirpath)
                % close

                % energy comparison for both configurations %<<<2
                figure
                x = t_pulsestart;
                y = E1.*sw1 + E2.*not(sw1);
                y2 = E2.*sw1 + E1.*not(sw1);
                plot(x, y, '-r', x, y2, '-b', x, y2-y, '-k' )
                legend('E_1', 'E_2', 'E_2-E_1')
                title([num2str(groupindex) ' - Energy of pulses, both currents, both configurations'])
                xlabel('time of pulse start (s)')
                ylabel('Energy (J)')
                saveplot(sprintf('%05d-energy_config_1_2', groupindex), dirpath)
                close

                % resistance %<<<2
                figure
                hold on
                plot(t_pulsestart, Ra(1,:), '-y', t_pulsestart, sgolayfilt(Ra(1,:), 2, 21), '-y', 'linewidth',2)
                plot(t_pulsestart, Ra(2,:), '-r', t_pulsestart, sgolayfilt(Ra(2,:), 2, 21), '-r', 'linewidth',2)
                plot(t_pulsestart, Rb(1,:), '-g', t_pulsestart, sgolayfilt(Rb(1,:), 2, 21), '-g', 'linewidth',2)
                plot(t_pulsestart, Rb(2,:), '-b', t_pulsestart, sgolayfilt(Rb(2,:), 2, 21), '-b', 'linewidth',2)
                hold off
                legend( 'Ra1', 'filt(Ra1)', 'Ra2', 'filt(Ra2)', 'Rb1', 'filt(Rb1)', 'Rb2', 'filt(Rb2)')
                title([num2str(groupindex) ' - Resistance of braking rheostats'])
                xlabel('time of pulse start (s)')
                ylabel('R (Ohm)')
                saveplot(sprintf('%05d-rheostats_1_2', groupindex), dirpath)
                close

        end
end % plots
