function [H_exp, H_an, LBF_exp, LBF_an, T_0] = P1_SSTD(expData)
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

%% Define Material Properties
    b = strsplit(expData.name,'_'); % gives a cell array (b) that is 1x3
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

%% Analytical
H_an = (expData.volts * expData.amps * (1/1000)) / (area_Rod * k_mat);

    H_exp = 0;
    LBF_exp = 0;
    LBF_an = 0;
    T_0 = 0;
end