% Extracts preselected part of signals so one can work with shorter part and save memory
%
% Based on readme.txt in directory "Alessandria-Novara-01-25-19-7.40pm"
% Vupf= Upstream filter Voltage (Pantograph Voltage)
% Vdwnf = Downstream filter Voltage
% Vhf = Voltage between the rheostat at lowest voltage  and ground
% Ip = Pantograph current
% HopA = Current in highest rheostat measured with an HOP800 transducer
% HopB = Current in lowest rheostat measured with an HOP800 transducer
% IrogA = Current in highest rheostat measured with a Rogowski coil transducer
% IrogB = Current in lowest rheostat measured with a Rogowski coil transducer

ids = 63000000;
ide = 67000000;

varnms = {'Vupf', 'Vdwnf', 'HopA', 'HopB', 'IrogA', 'IrogB', 'Vhf', 'Vupf'}
for i = 1:length(varnms)
        load([varnms{i} '.mat'])
        eval([varnms{i} ' = ' varnms{i} '(' num2str(ids) ':' num2str(ide) ');']);
        save('-v7', [varnms{i} '-selection.mat'], varnms{i});
        eval(['clear ' varnms{i}]);
endfor

% ids = 80039700;
% ide = 82437300;
% 
% load('Vupf.mat');
% Vupf = Vupf(ids:ide);
% load('Vdwnf.mat');
% Vdwnf = Vdwnf(ids:ide);
% 
% load('HopA.mat');
% HopA = HopA(ids:ide);
% load('HopB.mat');
% HopB = HopB(ids:ide);
% load('Ip.mat');
% Ip = Ip(ids:ide);
% load('IrogA.mat');
% IrogA = IrogA(ids:ide);
% load('IrogB.mat');
% IrogB = IrogB(ids:ide);
% load('Vdwnf.mat');
% Vdwnf = Vdwnf(ids:ide);
% load('Vhf.mat');
% Vhf = Vhf(ids:ide);
% load('Vupf.mat');
% Vupf = Vupf(ids:ide);
% 
% % remove undwanted variables:
% clear ans
% 
% save('-v4', 'example.mat')
