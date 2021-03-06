%% === Calculate energy of a braking group ====================
% A middle between pulses of first current waveform is found. At these points the current waveforms
% are split. The energy is calculated as integral in a whole part between these two split points.
% The two voltages are switched at each split for two currents. For the first pulse the first
% voltage is applied to first current -> configuration 1. For the first pulse the second voltage is
% applied to first current -> configuration 2. 

function [E] = energy(groupindex, fs, triglvl, dirpath, plots);

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
id_pulsestart = find(posslopes);

% get indexes of negative slopes (rise down):
negslopes = slopes;
negslopes(negslopes == 1) = 0;
id_pulseend = find(negslopes);

% get lengths between all slopes (in samples): 
lengths = diff(idslopes);

% get lengths between positive and negative slopes - when voltage is high (in samples): 
tmp = ismember(idslopes, id_pulsestart);
lengthson = lengths(tmp(1:end-1));

% get lengths between negative and positive slopes - when voltage is low (in samples): 
tmp = ismember(idslopes, id_pulseend);
lengthsoff = lengths(tmp(1:end-1));

% get lengths between positive slopes - periods of pulses (in samples): 
periods = diff(id_pulsestart);

% display various info
disp(['Length of signal: ' num2str(length(Ia)) ' samples, duration: ' num2str(length(Ia)./fs) ' s.'])
if isempty(periods)
        disp('No detected pulses. Only noise data');
        id_pulsestart = 1;
        id_pulsestartO = 1;
        splitidx1 = 1;
        splitidx2 = length(Ia);
else
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

        %% --- Get safety negative time offset of pulse starts -------------------- %<<<2
        offset = fix(min(lengthsoff)./2);
        disp(['Minimal off length is ' num2str(min(lengthsoff)) ' samples. All pulsestarts will be moved back by ' num2str(offset) ' samples.'])
        id_pulsestartO = id_pulsestart - offset;
        % fix first negative index:
        if id_pulsestartO(1) < 1
                id_pulsestartO(1) = 1;
        end
        % XXX now here could be some overlap, if id_pulsestart = [1 20 120 ] and offset = 50,
        % than one gets id_pulsestartO = [1 -30 70 ...] - how to solve?
        t_pulsestart = t(id_pulsestart);
        t_pulsestartO = t(id_pulsestartO);

        %% --- Fit time starts by line-------------------- %<<<2
        [P, S] = polyfit ([1:length(t_pulsestart)], t_pulsestart, 1);
        % pulse start varies easily by more than width
        % at the end some pulses are missing

        % indexes where split will happen:
        splitidx1 = id_pulsestartO; % index of start of split
        splitidx2 = splitidx1;             % index of end of split
        % add index for start of waveform:
        splitidx1 = [1 splitidx1];
        % add index for end of waveform:
        splitidx2(end+1) = length(Ia);
end

%% --- Split signal into individual current and nocurrent pulses --------------------
disp('Splitting signal...')

% make split:
Ia_splits = arrayfun( @(r1,r2) Ia(r1:r2), splitidx1, splitidx2, 'uni', false );
Ib_splits = arrayfun( @(r1,r2) Ib(r1:r2), splitidx1, splitidx2, 'uni', false );
Vhf_splits = arrayfun( @(r1,r2) Vhf(r1:r2), splitidx1, splitidx2, 'uni', false );
Vdsfhf_splits = arrayfun( @(r1,r2) Vdsfhf(r1:r2), splitidx1, splitidx2, 'uni', false );

%% --- Calculate power and energy -------------------- %<<<1
if isempty(periods)
        sw1 = 1;
else
        % prepare indexing vectors for two powers 
        % i.e. make series 1,0,1,0,... (based on square.m from statistics package):
        sw1 = [1:length(splitidx1)];
        sw1 = 0.5*sw1;
        sw1 = (sw1 - floor(sw1)>=0.5);
end

% multiply current and voltage for configuration 1 (a*hf, b*dsfhf)

%                pulse 1  pulse 2     pulse 3
%         P_1:   Ia*Vhf               Ia*Vhf
%                Ib*Vdsf-Vhf          Ib*Vdsf-Vhf
%                       Ia*Vdsf-Vhf
%                       Ib*Vhf
%
%         P_2:   Ia*Vdsf-Vhf
%                Ib*Vhf
%                    +
%                       Ia*Vhf
%                       Ib*Vdsf-Vhf
Phf_a    = Ia.*Vhf;
Phf_b    = Ib.*Vhf;
Pdsfhf_a = Ia.*Vdsfhf;
Pdsfhf_b = Ib.*Vdsfhf;
% split power waveforms
Phf_a_splits    = arrayfun( @(r1,r2) Phf_a(r1:r2),    splitidx1, splitidx2, 'uni', false );
Phf_b_splits    = arrayfun( @(r1,r2) Phf_b(r1:r2),    splitidx1, splitidx2, 'uni', false );
Pdsfhf_a_splits = arrayfun( @(r1,r2) Pdsfhf_a(r1:r2), splitidx1, splitidx2, 'uni', false );
Pdsfhf_b_splits = arrayfun( @(r1,r2) Pdsfhf_b(r1:r2), splitidx1, splitidx2, 'uni', false );
% calculate energy in every pulse:
Ehf_a    = 1./fs.*cellfun(@trapz, Phf_a_splits,    'uni', true);
Ehf_b    = 1./fs.*cellfun(@trapz, Phf_b_splits,    'uni', true);
Edsfhf_a = 1./fs.*cellfun(@trapz, Pdsfhf_a_splits, 'uni', true);
Edsfhf_b = 1./fs.*cellfun(@trapz, Pdsfhf_b_splits, 'uni', true);
% calculate total energy for configuration 1 (a*hf, b*dsfhf):
E(1,1) = sum(Ehf_a.*sw1) + sum(Edsfhf_b.*sw1) + sum(Edsfhf_a.*not(sw1)) + sum(Ehf_b.*not(sw1));
E(2,1) = sum(Edsfhf_a.*sw1) + sum(Ehf_b.*sw1) + sum(Ehf_a.*not(sw1)) + sum(Edsfhf_b.*not(sw1));

disp(['Error between two energies (E2-E1)/E2 (%): ' num2str((E(2)-E(1))/E(1).*100)]);
disp(['Error between two energies (E1-E2)/E1 (%): ' num2str((E(1)-E(2))/E(2).*100)]);

%% --- 2DO -------------------- XXX %<<<1
% calculate average offset level?

%% --- Plotting -------------------- %<<<1
if plots
        % plot current with lines showing splitting - this figure is challenging the hardware, use only if
        % needed!
        figure
        hold on
        if ~isempty(periods)
                x = repmat(t_pulsestartO, 2, 1);
                yl = [min(Ia) max(Ia)];
                y = repmat(yl', 1, length(t_pulsestartO));
                for i = 1:length(t_pulsestartO)
                        plot([t_pulsestartO(i) t_pulsestartO(i)], yl, '-r')
                end
                title(sprintf('Gr. %d - Current waveform Ia\nmean: %g A, std: %g A\nred: pulse start, green: pulse end', groupindex, mean(Ia), std(Ia)))
        else
                title(sprintf('Gr. %d - Current waveform Ia\nmean: %g A, std: %g A', groupindex, mean(Ia), std(Ia)))
        endif
        plot(t, Ia, '-b')
        hold off
        xlabel('time (s)')
        ylabel('I (A)')
        saveplot(sprintf('%05d-current_Ia', groupindex), dirpath)
        close

        if ~isempty(periods)
                % plot selected current pulses, so user can visually check if splitting occured correctly:
                figure
                hold on
                noofpulses = ceil(length(Ia_splits)./100);
                lens = [];
                for i = fix(linspace(1, length(Ia_splits), noofpulses))
                        plot(Ia_splits{i})
                        lens(end+1) = length(Ia_splits{i});
                end
                hold off
                title(sprintf('Gr. %d - Selected current pulses, Ia', groupindex))
                xlabel('time (samples), zero is set to pulse start')
                ylabel('current (A)')
                xlim([1 min(lens)])
                saveplot(sprintf('%05d-selected_current_pulses_Ia', groupindex), dirpath)
                close

                % errors of pulse starts from ideal line
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

                % energy in splits - pulses
                figure
                x = [1 t_pulsestart];
                y = Ehf_a.*sw1;
                x(y == 0) = [];
                y(y == 0) = [];
                x2 = [1 t_pulsestart];
                y2 = Edsfhf_a.*not(sw1);
                x2(y2 == 0) = [];
                y2(y2 == 0) = [];
                plot(x, y, 'xr', x2, y2, 'xb')
                legend('Ehf\_a', 'Edsfhf\_a')
                title([num2str(groupindex) ' - Energy of pulses, current a, configuration 1'])
                xlabel('time of pulse start (s)')
                ylabel('Energy (J)')
                saveplot(sprintf('%05d-energy_a_config_1', groupindex), dirpath)
                close

                % energy comparison for both configurations
                figure
                x = [1 t_pulsestart];
                y = Ehf_a.*sw1 + Edsfhf_a.*not(sw1);
                y2 = Edsfhf_a.*sw1 + Ehf_a.*not(sw1);
                plot(x, y, '-r', x, y2, '-b', x, y2-y, '-k' )
                legend('E_1', 'E_2', 'E_2-E_1')
                title([num2str(groupindex) ' - Energy of pulses, current a, both configurations'])
                xlabel('time of pulse start (s)')
                ylabel('Energy (J)')
                saveplot(sprintf('%05d-energy_a_config_1_2', groupindex), dirpath)
                close

        end
end % plots
