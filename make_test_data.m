% script generate testing data for calculate.m
% the signal contain:
% - every braking pulse 5 samples long
% - between pulses 20 samples
% - every group 15 pulses - intentionally not divisible by 2
% - at start, between groups and at end 3 seconds of signal
% - three groups of braking pulses
%
% variables IrogA, IrogB, Vdwnf, Vhf
% total energy:
% E_1 
% E_2 

% sampling frequency (Hz):
fs = 37.5e3;

Ia_lvl = 1000; % A
Ib_lvl = 1100; % A
Vdsf_lvl = 4000; % V
Vhf_lvl = 1800; % V

% braking pulse:
pulse = [ones(1, 5) zeros(1, 20)];
% braking group:
group = repmat(pulse, 1, 15);
% delay between groups:
delay = zeros(1, round(3*fs));
% start/end delay:
sedelay = zeros(1, round(0.1*fs));
% current variables:
Ia = Ia_lvl.*[delay group delay group delay group delay];
Ib = Ib_lvl.*[delay group delay group delay group delay];
Vdsf = Vdsf_lvl.*ones(size(Ia));
Vhf = Vhf_lvl.*ones(size(Ia));

% rename, transpose and save to keep measurement data format:
IrogA = Ia;
save('-v7', fullfile('testdata', 'IrogA.mat'), 'IrogA');
IrogB = Ib';
save('-v7', fullfile('testdata', 'IrogB.mat'), 'IrogB');
Vdwnf = Vdsf;
save('-v7', fullfile('testdata', 'Vdwnf.mat'), 'Vdwnf');
save('-v7', fullfile('testdata', 'Vhf.mat'), 'Vhf');

% calculation of energy:
% no of pulses  x  current  x  voltage  x  no of samples  /  by sampling frequency
Eg_1 = 3*(Ia_lvl.*Vhf_lvl.*5./fs) + 2*(Ib_lvl.*(Vdsf_lvl-Vhf_lvl).*5./fs)
Eg_2 = 2*(Ia_lvl.*Vhf_lvl.*5./fs) + 3*(Ib_lvl.*(Vdsf_lvl-Vhf_lvl).*5./fs)
E_1 = 2*Eg_1 + Eg_2
E_2 = 2*Eg_2 + Eg_1
disp(['difference 1: ' num2str((E_2 - E_1)./E_1*100) ' %']);
disp(['difference 2: ' num2str((E_1 - E_2)./E_2*100) ' %']);

figure; hold on; plot(Ia); plot(Ib); legend('Ia','Ib'); hold off
