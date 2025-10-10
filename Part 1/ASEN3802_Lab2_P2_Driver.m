clear all;
clc;
close all;
% Heat Conduction Lab
% Contributor(s): 
% Jackson Hughes
% Grant Bechtel
% Course: ASEN3802
% Date: 10/5/2025


data_Files = dir('*mA');
for i=1:length(data_Files)
    b = strsplit(data_Files(i).name,'_'); % gives a cell array (b) that is 1x3
    % {'material','voltsV','ampsmA'} -- now split by 'V' and 'mA'
    v = strsplit(b{2},'V'); % volts are always in the second portion
    ampval = strsplit(b{3},'mA'); % amps are always in the third portion
    volts(i) = str2num(v{1}); % convert string to number (vector)
    amps(i) = str2num(ampval{1});
    content_Files = readmatrix(data_Files(i).name);

    % read data into a struct
    expData(i).name = data_Files(i).name;
    expData(i).volts = volts(i);
    expData(i).amps = amps(i);
    expData(i).values = content_Files;
    clear content_Files;
end
clear b ampval v i data_Files volts amps content_Files;
%% Material Properties
in_to_m = 0.0254;

rho_Aluminum = 2810; % kg/m^3
rho_Brass = 8500; 
rho_Steel = 8000; 

cp_Aluminum = 960; % J/(kg*K)
cp_Brass = 380; 
cp_Steel = 500; 

k_Aluminum = 130; % W/(m*K)
k_Brass = 115;
k_Steel = 16.2;
%% task1
