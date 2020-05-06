%% === Calculate energy of a braking sequence ====================
function [E] = energy(groupindex, fs, triglvl, dirpath, plots);

varnms = {'Ia', 'Ib', 'Vdsf', 'Vhf'};
% generate all filenames
for j = 1:length(varnms)
        fng{j} = fullfile(dirpath, [varnms{j} sprintf('-%05d', groupindex) '.mat']);
        load(fng{j});
end

%% --- Other required quantities --------------------
% time axis:
t = [1:length(Ia)]./fs - 1/fs;
% half-filter voltage:
% complementary half-filter voltage:
Vdsfhf = Vdsf - Vhf;

%% --- Finding indexes of pulse starts and pulse ends and pulse lengths --------------------
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

        %% --- Get safety negative time offset of pulse starts --------------------
        offset = fix(min(lengthsoff)./2);
        disp([sprintf('\n') 'Minimal off length is ' num2str(min(lengthsoff)) ' samples, so all pulsestarts will be moved back by ' num2str(offset) ' samples.'])
        id_pulsestartO = id_pulsestart - offset;
        % fix first negative index:
        if id_pulsestartO(1) < 1
                id_pulsestartO(1) = 1;
        end
        % XXX now here could be some overlap, if id_pulsestart = [1 20 120 ] and offset = 50,
        % than one gets id_pulsestartO = [1 -30 70 ...] - how to solve?
        t_pulsestart = t(id_pulsestart);
        t_pulsestartO = t(id_pulsestartO);

        %% --- Fit time starts by line--------------------
        [P, S] = polyfit ([1:length(t_pulsestart)], t_pulsestart, 1);
        % pulse start varies easily by more than width
        % at the end some pulses are missing

        %% --- Split signal into individual current pulses --------------------
        disp('Splitting signal...')
        % indexes where split will happen:
        splitidx1 = id_pulsestartO; % index of start of split
        splitidx2 = splitidx1;             % index of end of split
        % add index for start of waveform:
        splitidx1 = [1 splitidx1];
        % add index for end of waveform:
        splitidx2(end+1) = length(Ia);
end

% make split:
Ia_splits = arrayfun( @(r1,r2) Ia(r1:r2), splitidx1, splitidx2, 'uni', false );
Ib_splits = arrayfun( @(r1,r2) Ib(r1:r2), splitidx1, splitidx2, 'uni', false );
Vhf_splits = arrayfun( @(r1,r2) Vhf(r1:r2), splitidx1, splitidx2, 'uni', false );
Vdsfhf_splits = arrayfun( @(r1,r2) Vdsfhf(r1:r2), splitidx1, splitidx2, 'uni', false );

%% --- Calculate power and energy --------------------
disp('Power and energy...')
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

%% --- 2DO -------------------- XXX
% calculate average offset level?

%% --- Plotting --------------------
if plots
        disp('Plotting...')
        % % plot current with lines showing splitting - this is figure challenging the hardware, use only if
        % % needed!
        % if ~isempty(periods)
        %         figure
        %         hold on
        %         x = repmat(t_pulsestartO, 2, 1);
        %         yl = [min(Ia) max(Ia)];
        %         y = repmat(yl', 1, length(t_pulsestartO));
        %         for i = 1:length(t_pulsestartO)
        %                 plot([t_pulsestartO(i) t_pulsestartO(i)], yl, '-r')
        %         end
        %         plot(t, Ia, '-b')
        %         hold off
        %         title('Current waveform')
        %         xlabel('time (s)')
        %         ylabel('I (A)')
        %         saveplot(sprintf('%05d-current_Ia', groupindex), dirpath)
        %         close
        % endif

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
        title([num2str(groupindex) ' - Selected current pulses, Ia'])
        xlabel('time (samples), zero is arbitrary')
        ylabel('current (A)')
        xlim([1 min(lens)])
        saveplot(sprintf('%05d-selected_current_pulses_Ia', groupindex), dirpath)
        close

        % % 4 selected consecutive pulses with strange behaviour
        % % probably not needed anymore XXX
        % figure
        % plot(t_pulsestart(1145:1155), Ia_pulse_int(1145:1155))
        % ylabel('I*t (A*s)')
        % xlabel('time (s)')
        % title('Integral of selected current pulses')
        % 
        % figure
        % hold on
        % plot(Ia_splits{1148})
        % plot(Ia_splits{1149})
        % plot(Ia_splits{1150})
        % plot(Ia_splits{1151})
        % hold off
        % title('Selected current pulses - 4 consecutive. Strange behaviour')
        % legend(['pulse 1148 at time ' num2str(t_pulsestart(1148))], ...
        %        ['pulse 1149 at time ' num2str(t_pulsestart(1149))], ...
        %        ['pulse 1150 at time ' num2str(t_pulsestart(1150))], ...
        %        ['pulse 1151 at time ' num2str(t_pulsestart(1151))])
        % xlabel('time (samples), zero is arbitrary')
        % xlabel('current (A)')
        % close

        if ~isempty(periods)
                % % time of pulse starts fitted by line
                % figure
                % hold on
                % plot(1:length(t_pulsestart), t_pulsestart, '-')
                % plot(1:length(t_pulsestart), polyval(P, 1:length(t_pulsestart)), 'x')
                % xlabel('index')
                % ylabel('time of pulse start (s)')
                % legend('time of pulse start', 'fit line')
                % title('fitting of pulse start times')
                % saveplot(sprintf('%05d-pulse_start_times', groupindex), dirpath)

                % errors of pulse starts from ideal line
                % reveals gaps between pulses
                figure
                hold on
                plot(polyval(P, [1:length(t_pulsestart)]) - t_pulsestart, '-b')
                plot(xlim, 0 + 0.5*[min(lengthsoff) min(lengthsoff)]./fs, 'r--');
                plot(xlim, 0 - 0.5*[min(lengthsoff) min(lengthsoff)]./fs, 'r--');
                legend('error', '0 - 1/2 of minimal "off" length', '0 - 1/2 of minimal "off" length')
                title([num2str(groupindex) ' - Errors of pulse starts time from ideal line fit'])
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
                legend('Ehf_a', 'Edsfhf_a')
                title([num2str(groupindex) ' - Energy of pulses, current a, for configuration 1'])
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
                title([num2str(groupindex) ' - Energy of pulses, current a'])
                xlabel('time of pulse start (s)')
                ylabel('Energy (J)')
                saveplot(sprintf('%05d-energy_a_config_1_2', groupindex), dirpath)
                close

        end
end % plots

%% --- Old unused code, stored for reference --------------------
% plot(x./3600, (DI.y.v - vnom).*1e6, '-');
% title([manufacturer ' ' zrtype ', s.n. ' sn ', nominal voltage ' tap_str]);
% ax = gca;
% pos = get(ax, 'position');
% set(ax, 'position', pos.*[1 1.5 0.9 0.85]);
% yl = get(ax,'ylim');
% yt = get(ax,'ytick');
% h0 = get(ax,'children');
% hold on
% [ax,h1,h2] = plotyy(ax,0,0,tempx./3600,tempy./1e3);
% delete(h1)
% set(ax(1),'ycolor',get(h0,'color'),'ylim',yl,'ytick',yt)
% set(h2, 'linestyle','--');
% ylabel(ax(2), 'Resistance of internal temperature sensor (k\Omega)')
% legend('location','southeast');
% legend('Voltage', 'Resistance')
% legend('boxoff');
% xlabel('Time (hour)')
% ylabel('Deviation from nominal value (\mu V)')
% 
% %% --- Switching times --------------------
% % First switch time is based on finding middle of first peak ramp up and moving back in time a little
% % This is not good method, it is not noise resistant, but I have no idea what noise can happen so I
% % do not know what to filter or what to expect.
% mid = (max(Ia)-min(Ia))/2;
% idx = find(Ia > mid);
% % time of first peak ramp up middle:
% tsw0 = t(idx(1));
% % % move time of first peak back by half chopping period:
% tsw0 = tsw0 - 1/fch/2;
% 
% % suppose all peaks are higher than triglvl!!! XXX
% Ia = [...
%         0.*ones(1, 5)...
%         1.*ones(1, 2)...
%         0.*ones(1, 5)...
%         1.*ones(1, 2)...
%         0.*ones(1, 5)...
%         1.*ones(1, 2)...
%         0.*ones(1, 5)...
%         1.*ones(1, 2)...
%         0.*ones(1, 5)...
%         1.*ones(1, 2)...
% ]
% 
% %% --- Switching times --------------------
% % Finds out at which time switch between rheostat a-b occured based on the chopping frequency.
% 
% 
%         % % % get times of all subsequent switches:
%         % % num = (max(t) - tsw0).*fs;
%         % % tsw = [1:num]./fch + tsw0;
%         % % % get indexes at chopping times
%         % % % nearest is for getting rounded values of indexes
%         % % idsw = interp1(t, [1:length(t)], tsw, 'nearest');
% 
% % make square function (0/1)
% %   based on square.m from package statistics of GNU Octave multiply by 0.5 because period should be
% %   2 times higher than chopping frequency so half period of square function covers one peak of
% %   current
% tmp = 0.5*fch*(t - tsw0);
% S = ones(size(tmp));
% S(tmp-floor(tmp) >= 0.5) = 0;
% S = max(Ia).*S;
% 
% %% --- Check proper sections --------------------
% 
% 
% 
% 
% id = 1;
% idprev = 1;
% finished = 0;
% onpeak = 0;
% do
%         % find next peak start
%         id = find(Ia(idprev, end));
%         id = id + idprev;
%         if isempty(id)
%                 finished = 1; 
%         end
%         % check slope
%         slope = (Ia(id+1) - Ia(id)) > 0;
%         if (onpeak & slope) | (~onpeak & ~slope)
%                 error(['onpeak and positive slope or not on peak and negative slope! id:' num2str(id)])
%         end
%         % measure distance from last peak
%         dist(end+1) = id - idprev;
%         if onpeak
%                 % calculate cumulative sum of current peak:
%                 intIp(end+1) = cumsum(Ia(idprev, id));
%                 % calculate power of the peak using a voltage
%                 Phf(end+1) = Ia(idprev, id).*Vhf(idprev, id);
%                 % calculate power of the peak using b voltage
%                 Pdsfhf(end+1) = Ia(idprev, id).*Vdsfhf(idprev, id);
%         % prepare for next iteration:
%         idprev = id;
% until finished
% 
% 
% 
% % CASY PREPNUTI se rozchazeji, vzorkovaci frce a fch nemaji stejnou casovou zakladnu. co s tim? kazdy
% % jednotlivy prechod zjistit samostatne! jenze co kdyz je vypadek v chopovani? spis je potreba
% % prepocitat casove zakladny!
% 
% %% --- Plot current with switching lines --------------------
% % XXX docasne, nemam vsechny veliciny
% figure
% hold on
% plot(t, Ia, 'r')
% plot(t,S,'b')
% % % yl = ylim;
% % % for i = 1:20
% % %         plot([tsw(i) tsw(i)], yl, '-k')
% % % end
% hold off
% xlim([4.72 4.9])

end % function
