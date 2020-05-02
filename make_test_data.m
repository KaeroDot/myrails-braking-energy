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
%
% Current (Ia and Ib):
%
%    3 seconds         15 pulses                3 seconds         15 pulses                3 seconds         15 pulses                3 seconds   
% <-------------><--------------------------><-------------><--------------------------><-------------><--------------------------><------------->
%                1      2             15                    16     17            30                    31     32            45
%                 _      _             _                     _      _             _                     _      _             _                    
%                | |    | |           | |                   | |    | |           | |                   | |    | |           | |                   
% _______________| |____| |____......_| |___________________| |____| |____......_| |___________________| |____| |____......_| |___________________
%                <-><-->
%                5 samples on state
%                   20 samples off state
% Power:
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
%
% Energy: P*samples/fs
%

%% quantity values:
% sampling frequency (Hz):
fs = 37.5e3;
% currents:
Ia_lvl = 1000; % A
Ib_lvl = 1100; % A
% voltages:
Vdsf_lvl = 4000; % V
Vhf_lvl = 1800; % V

%% calculate energy:
% no of pulses  x  current  x  voltage  x  no of samples  /  by sampling frequency
% Eg_1 = 15*Ia_lvl.* Vhf_lvl          .*5./fs + 15*Ib_lvl.* Vhf_lvl          .*5./fs
Ea_1 =  8*Ia_lvl.* Vhf_lvl          .*5./fs +  7*Ia_lvl.*(Vdsf_lvl-Vhf_lvl).*5./fs;
Eb_1 =  8*Ib_lvl.*(Vdsf_lvl-Vhf_lvl).*5./fs +  7*Ib_lvl.*Vhf_lvl           .*5./fs;
Eg_1 = Ea_1 + Eb_1
Ea_2 =  8*Ia_lvl.*(Vdsf_lvl-Vhf_lvl).*5./fs +  7*Ia_lvl.*Vhf_lvl           .*5./fs;
Eb_2 =  8*Ib_lvl.* Vhf_lvl          .*5./fs +  7*Ib_lvl.*(Vdsf_lvl-Vhf_lvl).*5./fs;
Eg_2 = Ea_2 + Eb_2

% Eg_2 = 15*Ia_lvl.*(Vdsf_lvl-Vhf_lvl).*5./fs + 15*Ib_lvl.*(Vdsf_lvl-Vhf_lvl).*5./fs
E_1 = 3*Eg_1
E_2 = 3*Eg_2
disp(['difference 1: ' num2str((E_2 - E_1)./E_1*100) ' %']);
disp(['difference 2: ' num2str((E_1 - E_2)./E_2*100) ' %']);

%% signal generation
% braking pulse:
pulse = [ones(1, 5) zeros(1, 20)];
% braking group:
group = repmat(pulse, 1, 15);
% delay between groups:
delay = zeros(1, round(3*fs));
% whole signal:
sig = [delay group delay group delay group delay];
% current variables:
Ia = Ia_lvl.*sig;
Ib = Ib_lvl.*sig;
Vdsf = Vdsf_lvl.*ones(size(Ia));
Vhf = Vhf_lvl.*ones(size(Ia));

%% rename to keep measurement data format and save:
% check direcotory:
if exist('testdata') ~= 7
        mkdir('testdata')
endif
IrogA = Ia;
save('-v7', fullfile('testdata', 'IrogA.mat'), 'IrogA');
IrogB = Ib';
save('-v7', fullfile('testdata', 'IrogB.mat'), 'IrogB');
Vdwnf = Vdsf;
save('-v7', fullfile('testdata', 'Vdwnf.mat'), 'Vdwnf');
save('-v7', fullfile('testdata', 'Vhf.mat'), 'Vhf');

%% plot:
% figure; hold on; plot(Ia); plot(Ib); legend('Ia','Ib'); hold off
