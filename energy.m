%% === Calculate energy of a braking sequence ====================
function [E_1 E_2] = energy(groupindex, fs, fch, triglvl, dirpath, plotting);

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
        disp(["\nMinimal off length is " num2str(min(lengthsoff)) ' samples, so all pulsestarts will be moved back by ' num2str(offset) ' samples.'])
        id_pulsestartO = id_pulsestart - offset;

        t_pulsestart = t(id_pulsestart);
        t_pulsestartO = t(id_pulsestartO);

        %% --- Fit time starts by line--------------------
        [P, S, MU] = polyfit ([1:length(t_pulsestart)], t_pulsestart, 1);
        % pulse start varies easily by more than width
        % at the end some pulses are missing

        %% --- Split signal into individual current pulses --------------------
        disp('Splitting signal...')
        % indexes where split will happen:
        splitidx1 = id_pulsestartO(1:end); % index of start of split
        splitidx2 = splitidx1;             % index of end of split
        % add index for start of waveform:
        splitidx1 = [1 splitidx1];
        % add index for end of waveform:
        splitidx2(end+1) = length(Ia);
endif

% make split:
Ia_splits = arrayfun( @(r1,r2) Ia(r1:r2), splitidx1, splitidx2, 'uni', false );
Ib_splits = arrayfun( @(r1,r2) Ib(r1:r2), splitidx1, splitidx2, 'uni', false );
Vhf_splits = arrayfun( @(r1,r2) Vhf(r1:r2), splitidx1, splitidx2, 'uni', false );
Vdsfhf_splits = arrayfun( @(r1,r2) Vdsfhf(r1:r2), splitidx1, splitidx2, 'uni', false );

%% --- Calculate current integral --------------------
% XXX probably this section is not needed anymore
% calculate integral of currents of pulses (./fs is for proper scaling):
Ia_pulse_int = 1./fs.*cellfun(@trapz, Ia_splits, 'uni', true);
Ib_pulse_int = 1./fs.*cellfun(@trapz, Ib_splits, 'uni', true);

%% --- Calculate power and energy --------------------
disp('Power and energy...')
if isempty(periods)
        sw1 = 1;
else
        % prepare indexing vectors for two powers 
        % i.e. make series 1,0,1,0,... (based on square.m from statistics package):
        sw1 = [1:length(t_pulsestart)+1]; %XXX opravdu t_pulsestart tady?
        sw1 = 0.5*sw1;
        sw1 = (sw1 - floor(sw1)>=0.5);
endif

% multiply current and voltage for configuration 1 (a*hf, b*dsfhf)
Phf_1 = Ia.*Vhf;
Pdsfhf_1 = Ib.*Vdsfhf;
% split power waveforms
Phf_1_splits = arrayfun( @(r1,r2) Phf_1(r1:r2), splitidx1, splitidx2, 'uni', false );
Pdsfhf_1_splits = arrayfun( @(r1,r2) Phf_1(r1:r2), splitidx1, splitidx2, 'uni', false );
% calculate energy in every pulse:
Ehf_1 = 1./fs.*cellfun(@trapz, Phf_1_splits, 'uni', true);
Edsfhf_1 = 1./fs.*cellfun(@trapz, Pdsfhf_1_splits, 'uni', true);
% calculate total energy for configuration 1 (a*hf, b*dsfhf):
E_1 = sum(Ehf_1.*sw1) + sum(Edsfhf_1.*not(sw1))

% multiply current and voltage for configuration 2 (b*hf, a*dsfhf)
Phf_2 = Ib.*Vhf;
Pdsfhf_2 = Ia.*Vdsfhf;
% split power waveforms
Phf_2_splits = arrayfun( @(r1,r2) Phf_2(r1:r2), splitidx1, splitidx2, 'uni', false );
Pdsfhf_2_splits = arrayfun( @(r1,r2) Phf_2(r1:r2), splitidx1, splitidx2, 'uni', false );
% calculate energy in every pulse:
Ehf_2 = 1./fs.*cellfun(@trapz, Phf_2_splits, 'uni', true);
Edsfhf_2 = 1./fs.*cellfun(@trapz, Pdsfhf_2_splits, 'uni', true);
% calculate total energy for configuration 1 (b*hf, a*dsfhf):
E_2 = sum(Ehf_2.*sw1) + sum(Edsfhf_2.*not(sw1))

disp(['Error between two energies (E2-E1)/E2 (%): ' num2str((E_2-E_1)/E_1.*100)]);


%% --- 2DO --------------------
% calculate average offset level

%% --- Plotting --------------------
disp('Plotting...')

if plotting
        % % plot current with lines showing splitting - this is figure challenging the hardware, use only if
        % % needed!
        % figure
        % hold on
        % x = repmat(t_pulsestartO, 2, 1);
        % yl = [min(Ia) max(Ia)];
        % y = repmat(yl', 1, length(t_pulsestartO));
        % for i = 1:length(t_pulsestartO)
        %         plot([t_pulsestartO(i) t_pulsestartO(i)], yl, '-r')
        % endfor
        % plot(t_pulsestart, Ia, '-b')
        % hold off
        % title('Current waveform')
        % xlabel('time (s)')
        % ylabel('I (A)')

        % plot selected current pulses, so user can visually check if splitting occured correctly:
        figure
        hold on
        for i = fix(linspace(1, length(Ia_splits), 100))
                plot(Ia_splits{i})
        endfor
        hold off
        title('Selected current pulses, Ia')
        xlabel('time (samples), zero is arbitrary')
        xlabel('current (A)')

        % plot integral of current pulses
        % probably not needed anymore XXX
        figure
        plot(t, Ia);
        title('Integrals of current pulses')
        ax = gca;
        pos = get(ax, 'position');
        set(ax, 'position', pos.*[1 1.5 0.9 0.85]);
        yl = get(ax,'ylim');
        yt = get(ax,'ytick');
        h0 = get(ax,'children');
        hold on
        [ax,h1,h2] = plotyy(ax,0,0, t_pulsestart, Ia_pulse_int);
        delete(h1)
        set(ax(1),'ycolor',get(h0,'color'),'ylim',yl,'ytick',yt)
        % set(h2, 'linestyle','--');
        % set(h2, 'markerstyle','x');
        ylabel(ax(2), 'I*t (A*s)')
        xlabel('time (s)')
        ylabel('I (A)')

        % 4 selected consecutive pulses with strange behaviour
        % probably not needed anymore XXX
        figure
        plot(t_pulsestart(1145:1155), Ia_pulse_int(1145:1155))
        ylabel('I*t (A*s)')
        xlabel('time (s)')
        title('Integral of selected current pulses')

        figure
        hold on
        plot(Ia_splits{1148})
        plot(Ia_splits{1149})
        plot(Ia_splits{1150})
        plot(Ia_splits{1151})
        hold off
        title('Selected current pulses - 4 consecutive. Strange behaviour')
        legend(['pulse 1148 at time ' num2str(t_pulsestart(1148))], ...
               ['pulse 1149 at time ' num2str(t_pulsestart(1149))], ...
               ['pulse 1150 at time ' num2str(t_pulsestart(1150))], ...
               ['pulse 1151 at time ' num2str(t_pulsestart(1151))])
        xlabel('time (samples), zero is arbitrary')
        xlabel('current (A)')

        % time of pulse starts fitted by line
        figure
        hold on
        plot(1:length(t_pulsestart), t_pulsestart, '-')
        plot(1:length(t_pulsestart), S.yf, 'x')
        xlabel('index')
        ylabel('time of pulse start (s)')
        legend('time of pulse start', 'fit line')
        title('fitting of pulse start times')
        xlim([0 14])
        ylim([4.73 4.81])

        % errors of pulse starts from ideal line
        figure
        plot(S.yf - t_pulsestart)
        title('Errors of pulse starts time from ideal line fit')
        xlabel('pulse start (s)')
        ylabel('error from line fit (s)')

        % energy in splits - pulses
        figure
        plot(t_pulsestart, Ehf_1.*sw1, 'xr', t_pulsestart, Edsfhf_1.*not(sw1), 'xb')
        title('Energy of pulses for configuration 1')
        xlabel('time of pulse start (s)')
        ylabel('Energy (J)')
endif % plotting



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
% % suppose all peaks are higher than triglvl!!! XXXXXXXXX
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

