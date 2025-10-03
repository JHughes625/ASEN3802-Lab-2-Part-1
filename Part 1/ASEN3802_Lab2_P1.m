clear all;
close all;
clc;

data_Files=dir('*mA');

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
clear b ampval v i data_Files;

%% Material Properties
in_to_m = 0.0254;

rho_Aluminum = 2810; % kg/m^3
rho_Brass = 8500; 
rho_Steel = 8000; 

cp_Aluminium = 960; % J/(kg*K)
cp_Brass = 380; 
cp_Steel = 500; 

k_Aluminium = 130; % W/(m*K)
k_Brass = 115;
k_Steel = 16.2;

%% Experimental Data

x_0 = (1 + (3/8)) * in_to_m; % m
radius_Rod = 0.5 * in_to_m; % m
area_Rod = radius_Rod^2 * pi; % m^2
L_Rod = x_0 + (0.5 * in_to_m * 7) + (1 * in_to_m); % m 

T_0 = expData(1).values(1,2); % TEMP SOLN. find better value later
% Move this into a function later
% TC Positions:
TC_Data(1,:) = [x_0 ...
    (x_0 + 1*(0.5 * in_to_m)) ...
    (x_0 + 2*(0.5 * in_to_m)) ...
    (x_0 + 3*(0.5 * in_to_m)) ...
    (x_0 + 4*(0.5 * in_to_m)) ...
    (x_0 + 5*(0.5 * in_to_m)) ...
    (x_0 + 6*(0.5 * in_to_m)) ...
    (x_0 + 7*(0.5 * in_to_m)) ...
    ];
% Take steady state temp of each TC
for i = 2:length(expData(1).values(1,:))
    TC_Data(2,i-1) = max(expData(1).values(:,i));
end

[polycoeff_1, S_1] = polyfit(TC_Data(1,:),TC_Data(2,:),1);
exp_LBF_1 = polycoeff_1(1).*TC_Data(1,:) + polycoeff_1(2);

% H_an = Qdot / (A * k_material)

%% Analytical 

H_an(1)= (expData(1).volts * expData(1).amps * (1/1000)) / (area_Rod * k_Aluminium);
H_an(2)= (expData(2).volts * expData(2).amps * (1/1000)) / (area_Rod * k_Aluminium);
H_an(3)= (expData(3).volts * expData(3).amps * (1/1000)) / (area_Rod * k_Brass);
H_an(4)= (expData(4).volts * expData(4).amps * (1/1000)) / (area_Rod * k_Brass);
H_an(5)= (expData(5).volts * expData(5).amps * (1/1000)) / (area_Rod * k_Steel);

T_analytical_LBF(1,:) = H_an(1) .* TC_Data(1,:) + polycoeff_1(2);

%% Plots
figure
hold on;
plot(expData(1).values(:,1), expData(1).values(:,2));
plot(expData(1).values(:,1), expData(1).values(:,3));
plot(expData(1).values(:,1), expData(1).values(:,4));
plot(expData(1).values(:,1), expData(1).values(:,5));
plot(expData(1).values(:,1), expData(1).values(:,6));
plot(expData(1).values(:,1), expData(1).values(:,7));
plot(expData(1).values(:,1), expData(1).values(:,8));
plot(expData(1).values(:,1), expData(1).values(:,9));
xlabel('Time [s]');
ylabel('Temperature [C]')
title('Aluminum 25V 240mA Thermocouple');
legend('Thermocouple 1', 'Thermocouple 2', 'Thermocouple 3', 'Thermocouple 4', ...
    'Thermocouple 5', 'Thermocouple 6', 'Thermocouple 7', 'Thermocouple 8');
hold off;
figure
plot(TC_Data(1,:),exp_LBF_1(:));





