% Script calculates Energy Dissipated in Braking Rheostats in DC Railway Systems

% Quantity names based on abstract:
% Accurate Measurement of Energy Dissipated in Braking Rheostats in DC Railway Systems
% Helko van den Brom, Domenico Giordano, Danielle Gallo, Andreas Wank, Yljon Seferi
% hvdbrom@vsl.nl

%% CONFIGURATION %%%%%%%%%%%%%%%%%%%%%%%%
% sampling frequency (Hz):
fs = 37.5e3;
% directory with all data:
dirpath = 'exampledata';
dirpath = 'testdata';
% make plots? (0/1);
plots = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% --- Load data --------------------
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

%% --- Trigger level --------------------
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
% braking waveform, contains only 0, 1 %XXXXXXXXXXXx what about IrogB?
braking = Ia > triglvl;

%% --- Find braking groups --------------------
% A group of braking ends if no breaking for at least 1 second long.
% Reshape current into groups 1 seconds long and find indexes of groups when braking
% starts/stops.

% group duration in seconds:
group_duration = 1;
% number of samples in a group:
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

%% --- Plot braking groups --------------------
if plots
        figure
        hold on
        plot(groupsids, groups)
        yl = ylim;
        for i = 1:length(groups_start_id)
                plot([groups_start_id(i) groups_start_id(i)], yl)
        end
        hold off
        xlabel('Index of Ia')
        ylabel('Braking intensity (arbitrary)')
        title('Braking detection, braking groups, starts of groups')
        saveplot('braking_groups', dirpath);
        close
end

% clear quantities not needed anymore:
clear braking;

%% --- Divide data into groups and save data as new files --------------------
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

% load and splits for those variables and filenames that do not exist:
for j = 1:length(varnms)
        if ~all(isfile(j, :))
                % XXX do not load IrogA - already loaded
                load(fullfile(dirpath, [fnms{j} '.mat']));
                eval([fnms{j} 'all = ' fnms{j} '; clear ' fnms{j} ';']);
                for i = 1:length(groups_start_id)
                        if ~isfile(j, i)
                                % cuts variable to single group
                                % and rename it to be consistent with the paper
                                % and ensure row vectors:
                                eval([varnms{j} ' = ' fnms{j} 'all(' num2str(groups_start_id(i)) ':' num2str(groups_end_id(i)) ');']);
                                eval([varnms{j} ' = ' varnms{j} '(:)' char(39) ';']);
                                save('-v7', fng{j, i}, varnms{j});
                        end
                end
                eval(['clear ' varnms{j}]);
        end
end

%% --- Calculate energy for individual groups --------------------
E_1 = zeros(size(groups_start_id));
E_2 = zeros(size(groups_start_id));
for i = 1:length(groups_start_id)
        disp([' === group ' num2str(i) ' from ' num2str(length(groups_start_id)) ' ===']);
        % function [E_1 E_2] = energy(groupindex, fs, triglvl, dirpath, plotting);
        [E_1(i) E_2(i)] = energy(i, fs, triglvl, dirpath, plots);
end

disp([' === estimates for total data  ===']);
disp(['Total energy 1: ' num2str(sum(E_1)) ' J, total energy 2: ' num2str(sum(E_2)) ' J.']);
disp(['Error between two energies (E2-E1)/E2 (%):']);
(sum(E_2)-sum(E_1))./sum(E_1).*100
disp(['Error between two energies (E1-E2)/E1 (%):']);
(sum(E_1)-sum(E_2))./sum(E_2).*100
