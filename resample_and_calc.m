% script resample data and calculates energy for multiple resampling frequencies

tic
% sampling frequency of original data (Hz):
fs = 50e3;
% original data directory
origdatadir = './exampledata';
% resampled data directory prefix
dataprefix = '.';
% filenames to process:
% fnms = {'IrogA', 'IrogB', 'Vdwnf', 'Vhf'};
% varsold = {'IrogAvo', 'IrogBvo', 'Vdwnfvo', 'Vhfvo'};

% resampling coefficients:
% (last value is without resampling)
% P = [10:2:99, 100];
% Q = 100.*ones(size(P));
% % resampled frequencies:
% fsr = P./Q.*fs
Q = [1:20];
% Q = [1:4];
fsr = fs./Q

% common plot indexes/time:
pts = 792160./fs;
pte = 792300./fs;
% pts = 1008090./fs;
% pte = 1008490./fs;
% pts = 930800/fs;
% pte = 931100/fs;
% pts = 1234440/fs;
% pte = 1234760/fs;

% load data
load(fullfile(origdatadir, 'IrogA.mat'))
load(fullfile(origdatadir, 'IrogB.mat'))
load(fullfile(origdatadir, 'Vdwnf.mat'))
load(fullfile(origdatadir, 'Vhf.mat'))
% copy data
IrogAo = IrogA;
IrogBo = IrogB;
Vdwnfo = Vdwnf;
Vhfo = Vhf;

% loop over all resamples and prepare data
plotdata = {};
disp('preparing data...')
for i = 1:length(Q)
        % directory path
        dirpath{i} = fullfile(dataprefix, ['data_' num2str(i, '%04d')]);
        % prepare directory
        mkdir(dirpath{i});
        if i == 1
                % no resampling/decimating. just original data
                IrogA = IrogAo;
                IrogB = IrogBo;
                Vdwnf = Vdwnfo;
                Vhf = Vhfo;
        else
                % resample data
                % IrogA = resample(IrogAo, P(i), Q(i));
                % IrogB = resample(IrogBo, P(i), Q(i));
                % Vdwnf = resample(Vdwnfo, P(i), Q(i));
                % Vhf = resample(Vhfo, P(i), Q(i));
                IrogA = decimate(IrogAo, Q(i),1);
                IrogB = decimate(IrogBo, Q(i),1);
                Vdwnf = decimate(Vdwnfo, Q(i),1);
                Vhf =   decimate(Vhfo,   Q(i),1);
                % IrogA = IrogAo(1:Q(i):end); % this is equal to decimate(,,1)
                % IrogB = IrogBo(1:Q(i):end); % this is equal to decimate(,,1)
                % Vdwnf = Vdwnfo(1:Q(i):end); % this is equal to decimate(,,1)
                % Vhf =   Vhfo  (1:Q(i):end); % this is equal to decimate(,,1)
        endif
        % save data
        save('-v7', fullfile(dirpath{i}, 'IrogA.mat'), 'IrogA')
        save('-v7', fullfile(dirpath{i}, 'IrogB.mat'), 'IrogB')
        save('-v7', fullfile(dirpath{i}, 'Vdwnf.mat'), 'Vdwnf')
        save('-v7', fullfile(dirpath{i}, 'Vhf.mat'), 'Vhf')
        % get waveform for plotting:
        plotdata{i} = IrogA(fix(pts.*fsr(i)):ceil(pte.*fsr(i)));
end
        
% loop over all resamples and calculate data
for i = 1:length(dirpath)
        disp(['Energy estimation for fs = ' num2str(fsr(i))])
        % run calculation
        [tmp, tmp2] = calculate(dirpath{i}, fsr(i));
        % tmp is all possible combinations, so mean(tmp) is the mean estimate of energy E
        E(i) = mean(tmp);
        % similar for tmp2 -> uE
        uE(i) = max(tmp2);
endfor

% plot of energy and uncertainty:
figure
plot(fsr, E./E(1))
xlabel('resampled frequency (Hz)')
ylabel('energy, normalized')
title('Energy on sampling frequency')

figure
plot(fsr, uE./uE(1))
xlabel('resampled frequency (Hz)')
ylabel('uncertainty of energy, normalized')
title('Uncertinaty of energy on sampling frequency')

% plot of peaks:
figure
hold on
leg = {};
for i = 1:length(plotdata)
        plot([1:length(plotdata{i})]./fsr(i), plotdata{i})
        leg{i} = num2str(fsr(i));
end
xlabel('time (s)')
ylabel('voltage (V)')
legend(leg)
title('Selected peak for various sampling frequencies')

disp('total calculation time:')
toc
