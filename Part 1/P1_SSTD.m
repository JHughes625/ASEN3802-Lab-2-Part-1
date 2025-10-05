% Heat Conduction Lab
% Contributor(s): 
% Jackson Hughes
% Course: ASEN3802
% Date: 10/5/2025

function [H_exp, LBF_exp, H_an, LBF_an, T_0] = P1_SSTD(expData)
%% Material Properties
in_to_m = 0.0254;

rho_Aluminum = 2810; % [kg/m^3]
rho_Brass = 8500; 
rho_Steel = 8000; 

cp_Aluminum = 960; % [J/(kg*K)]
cp_Brass = 380; 
cp_Steel = 500; 

k_Aluminum = 130; % [W/(m*K)]
k_Brass = 115;
k_Steel = 16.2;

%% Experimental Data

x_0 = (1 + (3/8)) * in_to_m; % [m]
radius_Rod = 0.5 * in_to_m; % [m]
area_Rod = radius_Rod^2 * pi; % [m^2]
L_Rod = x_0 + (0.5 * in_to_m * 7) + (1 * in_to_m); % [m] 

%% Define Material Properties

% Take in name string and detect which material to use
    b = strsplit(expData.name,'_');
    if strcmp(b{1}, 'Aluminum') == 1
        rho_mat = rho_Aluminum;
        cp_mat = cp_Aluminum;
        k_mat = k_Aluminum;
    elseif strcmp(b{1}, 'Brass') == 1
        rho_mat = rho_Brass;
        cp_mat = cp_Brass;
        k_mat = k_Brass;
    elseif strcmp(b{1}, 'Steel') == 1
        rho_mat = rho_Steel;
        cp_mat = cp_Steel;
        k_mat = k_Steel;
    end
% Define locations of thermocouples
TC_Data(1,:) = [
    x_0 ...
    (x_0 + 1*(0.5 * in_to_m)) ...
    (x_0 + 2*(0.5 * in_to_m)) ...
    (x_0 + 3*(0.5 * in_to_m)) ...
    (x_0 + 4*(0.5 * in_to_m)) ...
    (x_0 + 5*(0.5 * in_to_m)) ...
    (x_0 + 6*(0.5 * in_to_m)) ...
    (x_0 + 7*(0.5 * in_to_m)) ...
    ];

%% Analytical

% assume T_0 is the same as the initial temp of the TC closest to cold end
T_0 = expData.values(1,2); % [C]

%use given formula to calc. H_an
H_an = (expData.volts * expData.amps * (1/1000)) / (area_Rod * k_mat); % [C/m]
LBF_an = H_an .* TC_Data(:) + T_0;

%% Experimental

for i = 2:length(expData(1).values(1,:))
    % Take the mean of the values in the range of the steady state length
    % of the data
    TC_Data(2,i-1) = mean(expData.values(150:length(expData.values)-2,i));
end
% use polyfit to generate a Line of Best Fit for the dataset
[polycoeff_1] = polyfit(TC_Data(1,:),TC_Data(2,:), 1);

% pass soln. values to function
H_exp = polycoeff_1(1); % [C/m]
LBF_exp = polycoeff_1(1).*TC_Data(1,:) + polycoeff_1(2);

end