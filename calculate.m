% Script calculates Energy Dissipated in Braking Rheostats in DC Railway Systems

% Quantity names based on abstract:
% Accurate Measurement of Energy Dissipated in Braking Rheostats in DC Railway Systems
% Helko van den Brom, Domenico Giordano, Danielle Gallo, Andreas Wank, Yljon Seferi
% hvdbrom@vsl.nl

clear all;

%% CONFIGURATION %%%%%%%%%%%%%%%%%%%%%%%%
% sampling frequency (Hz):
fs = 37.5e3;
% directory with all data:
dirpath = 'exampledata';
% dirpath = 'testdata';
% make plots? (0/1):
plots = 1;
% rewrite splitted data? (can take long time, has no sense if this script is not changed):
rewrite = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% --- Load data -------------------- %<<<1
% load example data
%       Based on readme.txt in directory "Alessandria-Novara-01-25-19-7.40pm"
%       Vupf= Upstream filter Voltage (Pantograph Voltage)
%       Vdwnf = Downstream filter Voltage
%       Vhf = Voltage between the rheostat at lowest voltage  and ground
%       Ip = Pantograph current
%       HopA = Current in highest rheostat measured with an HOP800 transducer
%       HopB = Current in lowest rheostat measured with an HOP800 transducer
%       IrogA = Current in highest rheostat measured with a Rogowski coil transducer
%       IrogB = Current in lowest rheostat measured with a Rogowski coil transducer
% loaded quantities are renamed to naming used in paper
load(fullfile(dirpath, 'IrogA.mat'));
Ia = IrogA(:)';
clear IrogA;

%% --- Trigger level -------------------- %<<<1
% get trigger level as a middle of two histogram maxima. Hopefully it is more noise resistant than
% taking half of maximum voltage.
% histogram:
[NN, XX] = hist(Ia, 10);
% get maximum:
[tmp max1id] = max(NN);
max1id = max1id(1);
% mask first maximum:
NN(max1id) = 0;
% get second maximum:
[tmp max2id] = max(NN);
max2id = max2id(1);
% trigger level as half distance between two most common voltage levels:
triglvl = abs(XX(max2id) - XX(max1id))./2;
% braking waveform, contains only 0, 1
% script suppose Ib has the same triglvl and the same position of pulses
braking = Ia > triglvl;

%% --- Find braking groups -------------------- %<<<1
% A group of braking ends if no breaking for at least 1 second long.
% Reshape current into groups 1 seconds long and find indexes of groups when braking
% starts/stops.

% group duration in seconds:
group_duration = 1;

group_samples = fix(group_duration.*fs);
% reshape braking into matrix, one column contain samples of one group:
len = fix(size(braking,2)./group_samples);
groups = reshape(braking(1:group_samples.*len), group_samples, len); 
% make sum through column (i.e. inside one group):
groups = sum(groups, 1);
% generate indexes of braking curve that denote starts of groups:
groupsids = [1:group_samples:length(braking)];
% add sum of the last part of the braking curve that do not contain full number of samples to make
groups(end+1) = sum(groups(groupsids(end):end));

% find braking starts by finding upward slopes:
slopes = diff(groups > 0);
slopes(slopes == -1) = 0;
groups_start = find(slopes);
% now the groups_start contains indexes of groups variable, we need to convert it to indexes of
% braking variable:
groups_start_id = groups_start.*group_samples - group_samples;
groups_end_id = [groups_start_id length(braking)];
groups_start_id = [1 groups_start_id];

%% --- Plot braking groups -------------------- %<<<1
if plots
        figure
        hold on
        plot(groupsids, groups)
        yl = ylim;
        xtmp =  repmat(groups_start_id, 2, 1);
        ytmp = repmat(yl', 1, size(groups_start_id, 2));
        plot(xtmp, ytmp, '-r');
        hold off
        legend('braking', 'group start')
        xlabel('time (samples)')
        ylabel('braking intensity based on Ia (arbitrary)')
        title('Braking detection, braking groups')
        saveplot('braking_groups', dirpath);
        close
end

% clear quantities not needed anymore:
clear braking;

%% --- Divide data into groups and save data as new files -------------------- %<<<1
% filenames to process:
fnms = {'IrogA', 'IrogB', 'Vdwnf', 'Vhf'}; % IrogA is already loaded as Ia
% variables to process (paper use different variable names than measured files)
varnms = {'Ia', 'Ib', 'Vdsf', 'Vhf'}; % IrogA is already loaded as Ia
% generate all filenames to check and write:
for j = 1:length(varnms)
        for i = 1:length(groups_start_id)
                fng{j, i} = fullfile(dirpath, [varnms{j} sprintf('-%05d', i) '.mat']);
                % check if file exist:
                isfile(j, i) = exist(fng{j, i});
        end
end

% rewrite of splitted data overridden?
if rewrite
        isfile = isfile.*0;
end

% load and splits for those variables and filenames that do not exist:
for j = 1:length(varnms)
        if ~all(isfile(j, :))
                if j == 1
                        % IrogA already loaded as Ia, just make renaming so rest of the loop can
                        % work in the same way as for other variables;
                        eval([fnms{j} ' = ' varnms{j} ';']);
                        eval(['clear ' varnms{j}]);
                else
                        % load variable from measured data:
                        load(fullfile(dirpath, [fnms{j} '.mat']));
                end
                % rename measured data to SOMETHINGall:
                eval([fnms{j} 'all = ' fnms{j} '; clear ' fnms{j} ';']);
                for i = 1:length(groups_start_id)
                        if ~isfile(j, i)
                                % cuts variable to single group
                                % rename it to be consistent with the paper
                                eval([varnms{j} ' = ' fnms{j} 'all(' num2str(groups_start_id(i)) ':' num2str(groups_end_id(i)) ');']);
                                % ensure row vectors:
                                eval([varnms{j} ' = ' varnms{j} '(:)' char(39) ';']);
                                % save variable into file:
                                save('-v7', fng{j, i}, varnms{j});
                        end
                end
                eval(['clear ' varnms{j}]);
                eval(['clear ' fnms{j} 'all']);
        end
end

%% --- Calculate energy for individual groups -------------------- %<<<1
E = zeros(2, size(groups_start_id, 2));
EN = zeros(2, size(groups_start_id, 2));
uE = zeros(2, size(groups_start_id, 2));
urE = zeros(2, size(groups_start_id, 2));
for i = 1:length(groups_start_id)
        disp([' === group ' num2str(i) ' from ' num2str(length(groups_start_id)) ' ===']);
        % [E(:,i)] = energy(i, fs, triglvl, dirpath, plots);
        [E(:,i) EPN(:,i) EN(:,i) urE(:,i)] = energy3(i, fs, triglvl, dirpath, plots);
        uE = urE.*E;
end

%% --- Calculate total energy -------------------- %<<<1
disp([' =================================']);
disp([' === estimates for total data  ===']);
% get all possibilities of E_1 and E_2
% (i.e. all possible combinations of setups in braking groups. We do not know which combination is
% correct one, so lets calculate mean and std)
% permutations with repetitions:
E = sum(permrep(E), 2);
disp(['Braking energy                                     : ' num2str(mean(E)) ' J.']);
% uncertainty components:
uE = sqrt(sum(permrep(uE).^2, 2));
disp(['Uncertainty caused by pulse noise fitting          : ' num2str(max(uE)) ' J (' num2str(max(uE)./mean(E).*100) ' %).']);
disp(['Uncertainty caused by unknown correct configuration: ' num2str(std(E)) ' J (' num2str(std(E)./mean(E).*100) ' %).']);
% do the same for energy of noise:
EPN = sum(permrep(EPN), 2);
EN = sum(permrep(EN), 2);
disp(['Noise during pulse energy: (' num2str(mean(EPN)) ' +- ' num2str(std(EPN)) ') J (' num2str(std(EPN)./mean(EPN).*100) ' %).']);
disp(['Ratio of noise during pulse energy and braking energy is ' num2str(mean(EPN)./mean(E).*100) ' %.']);
disp(['Total noise energy: (' num2str(mean(EN)) ' +- ' num2str(std(EN)) ') J (' num2str(std(EN)./mean(EN).*100) ' %).']);
disp(['Total noise energy is ' num2str(mean(EN)./(mean(E) + mean(EN))*100) ' % of total energy.']);
