clear all;
clc;
close all;

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
%clear b ampval v i data_Files;

for i=1:length(expData)
    [H_exp(i), H_an(i), LBF_exp(i), LBF_an(i), T_0(i)] = P1_SSTD(expData(i));
end

for j=1:length(expData) %TASK 2 Plotting
[g(j,:),M(j),x_L] = M_exp(expData(j).values);
figure()
plot(x_L,g(j,:),linewidth=1.3,color=[224/255, 115/255, 52/255])%plot exp. IC via LOBF
%plot(,,linewidth=1.3,color=[52/255, 144/255, 224/255]) %plot analytical initial conditions
strTitle = expData(i).name +" Initial Condition Comp.";
title(strTitle)
xlabel("Dist From Th1 (m)")
ylabel("Temp C")
legend("Exp. LOBF IC", "Analytical IC")
end
%M % M is here to print values to input into overleaf table.