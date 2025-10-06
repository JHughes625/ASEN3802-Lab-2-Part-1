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

%% Experimental Data

x_0 = (1 + (3/8)) * in_to_m; % m
radius_Rod = 0.5 * in_to_m; % m
area_Rod = radius_Rod^2 * pi; % m^2
L_Rod = x_0 + (0.5 * in_to_m * 7) + (1 * in_to_m); % m 
TC_Positions = [
    x_0 ...
    (x_0 + 1*(0.5 * in_to_m)) ...
    (x_0 + 2*(0.5 * in_to_m)) ...
    (x_0 + 3*(0.5 * in_to_m)) ...
    (x_0 + 4*(0.5 * in_to_m)) ...
    (x_0 + 5*(0.5 * in_to_m)) ...
    (x_0 + 6*(0.5 * in_to_m)) ...
    (x_0 + 7*(0.5 * in_to_m)) ...
    ];

for i=1:length(expData)
    [H_exp(i), LBF_exp(:,i), H_an(i), LBF_an(:,i), T_0(i)] = P1_SSTD(expData(i));
end



%% Plots & Tables

outputTable = table(T_0,H_an,H_exp)

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
hold on
plot(TC_Positions,LBF_an(:,1));
plot(TC_Positions,LBF_exp(:,1));
xlabel('Position [m]');
ylabel('Temperature [C]')
title('Aluminum 25V 240mA SSTD Comparison');
legend('[An]','[Exp]');
hold off

figure
hold on
plot(TC_Positions,LBF_an(:,2));
plot(TC_Positions,LBF_exp(:,2));
xlabel('Position [m]');
ylabel('Temperature [C]')
title('Aluminum  30V 290mA SSTD Comparison');
legend('[An]','[Exp]');
hold off

figure
hold on
plot(TC_Positions,LBF_an(:,3));
plot(TC_Positions,LBF_exp(:,3));
xlabel('Position [m]');
ylabel('Temperature [C]')
title('Brass  25V 237mA SSTD Comparison');
legend('[An]','[Exp]');
hold off

figure
hold on
plot(TC_Positions,LBF_an(:,4));
plot(TC_Positions,LBF_exp(:,4));
xlabel('Position [m]');
ylabel('Temperature [C]')
title('Brass  30V 285mA SSTD Comparison');
legend('[An]','[Exp]');
hold off

figure
hold on
plot(TC_Positions,LBF_an(:,5));
plot(TC_Positions,LBF_exp(:,5));
xlabel('Position [m]');
ylabel('Temperature [C]')
title('Steel 22V 203mA SSTD Comparison');
legend('[An]','[Exp]');
hold off

for j=1:length(expData) %TASK 2 Plotting
[g(j,:),M(j),x_L] = M_exp(expData(j).values);
figure()
hold on
plot(x_L,g(j,:),linewidth=1.3,color=[224/255, 115/255, 52/255])%plot exp. IC via LOBF
plot(TC_Positions - TC_Positions(1),expData(j).values(1,2:9),'kx') % overlay the data
yline(T_0(i),linewidth=1.3,color=[52/255, 144/255, 224/255])
ax = gca;  % get current axes handle
xlimits = ax.XLim; 
ylimits = ax.YLim;
xlim([xlimits(1)-0.01, xlimits(2)+0.01])
ylim([ylimits(1)-0.1, ylimits(2)+0.1])
strTitle = expData(j).name +" Initial Condition Comp.";
strTitleSave = "PlotIC_"+expData(j).name ;
title(strTitle, 'interpreter', 'none')
xlabel("Dist From Th1 (m)")
ylabel("Temp C")
legend("Exp. LOBF IC", "Exp. Data", "Analytical IC", location='best')
hold off
saveas(gcf,strTitleSave,'png')
end
%M % M is here to print values to input into overleaf table.