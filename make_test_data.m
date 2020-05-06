% script generate testing data for calculate.m
% the signal contain:
% - every braking pulse 5 samples long
% - between pulses 20 samples
% - every group 15 pulses - intentionally not divisible by 2
% - at start, between groups and at end 3 seconds of signal
% - three groups of braking pulses
% - additional offset current added to whole data
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
%                <-> 5 samples on state
%                   <--> 20 samples off state
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
% current offset (simulates nonzero current):
Io_lvl = 0.001; % A

disp(' === Testing data ===')
disp(['Current Ia: ' num2str(Ia_lvl) ' A, Ib: ' num2str(Ib_lvl) ' A.'])
disp(['Current offset Io: ' num2str(Io_lvl) ' A.'])
disp(['Voltage Vdsf: ' num2str(Vdsf_lvl) ' V, Vhf: ' num2str(Vhf_lvl) ' V.'])
disp('')

%% signal generation
% normalized braking pulse:
pulse = [ones(1, 5) zeros(1, 20)];
% braking group:
group = repmat(pulse, 1, 15);
% delay between groups:
delay = zeros(1, round(3*fs));
% whole normalized signal:
sig = [delay group delay group delay group delay];
disp(['Signal length: ' num2str(length(sig)) ', signal duration: ' num2str(length(sig)./fs) ' s.']);
% current waveforms with offset:
Ia = Ia_lvl.*sig + Io_lvl;
Ib = Ib_lvl.*sig + Io_lvl;
% voltage waveforms:
Vdsf  = Vdsf_lvl.*ones(size(sig));
Vhf   = Vhf_lvl.*ones(size(sig));

%% calculate energies:
% (index) is configuration 1 or configuration 2

% energy in pulses:
% no of pulses  x  current  x  voltage  x  no of samples  /  by sampling frequency
Ea(1) =  8*Ia_lvl.* Vhf_lvl          .*5./fs +  7*Ia_lvl.*(Vdsf_lvl-Vhf_lvl).*5./fs;
Ea(2) =  8*Ia_lvl.*(Vdsf_lvl-Vhf_lvl).*5./fs +  7*Ia_lvl.*Vhf_lvl           .*5./fs;
Eb(1) =  8*Ib_lvl.*(Vdsf_lvl-Vhf_lvl).*5./fs +  7*Ib_lvl.*Vhf_lvl           .*5./fs;
Eb(2) =  8*Ib_lvl.* Vhf_lvl          .*5./fs +  7*Ib_lvl.*(Vdsf_lvl-Vhf_lvl).*5./fs;
Eg = Ea + Eb;
disp('Energy of single group, only in pulses (without offset current):')
disp(['configuration 1: ' num2str(Eg(1)) ' J, configuration 2: ' num2str(Eg(2)) ' J.'])

% total pulse energy:
E = 3*Eg;
disp('Total energy, only in pulses (without offset current):')
disp(['configuration 1: ' num2str(E(1)) ' J, configuration 2: ' num2str(E(2)) ' J.'])
disp(['Difference between configurations (2-1)/1: ' num2str((E(2) - E(1))./E(1)*100) ' %']);
disp(['Difference between configurations (1-2)/2: ' num2str((E(1) - E(2))./E(2)*100) ' %']);

% energy in offset:
Eo(1) = Io_lvl.*Vhf_lvl.*length(sig)./fs;
Eo(2) = Io_lvl.*(Vdsf_lvl - Vhf_lvl)*length(sig)./fs;
% two times because two currents:
Eo = 2*Eo;
disp(['Offset level: ' num2str(Io_lvl) ' A, energy only in offset:']);
disp(['configuration 1: ' num2str(Eo(1)) ' J, configuration 2: ' num2str(Eo(2)) ' J.'])
disp('Ratio of energy in offset to energy in pulses:')
disp(['Offset conf. 1 / energy conf. 1: ' num2str(Eo(1)./E(1)*100) ' %, off. c. 2 / en. c. 2: ' num2str(Eo(2)./E(2)*100) ' %.']) 
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
