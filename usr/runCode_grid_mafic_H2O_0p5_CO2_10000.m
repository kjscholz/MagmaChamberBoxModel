clear all
close all
clc

addpath('../src/')

tic
%% Setting up initial parameters


%param.composition = 'silicic';
param.composition = 'mafic';

run set_properties
run set_numerical_params

save_output  = 1; % set to 1 to save output as .mat file
save_figures = 0; % set to 1 to save figures


%% create grid for chamber radius, volatile content, recharge rate, and depth

% chamber volume (km^3) -- note: enter the LOG of volume
range_log_volume_km3 = linspace(0,2,30);  % range of log volume in km3

% volatile content
range_water      = linspace(4,0.5,1); % h2o content (wt%)
range_co2        = linspace(100,10000,1); % co2 content (ppm)

% recharge rate (km3/yr) -- note: enter the LOG of recharge rate
range_log_vfr    = linspace(-4,-2,30); % log volume recharge rate (km3/yr)

% storage depth (m)
range_depth      = linspace(8e3,8e3,1); % chamber depth (m)


%% Run the code over the grid search

% do some conversions & create the grid
range_volume_km3 = 10.^range_log_volume_km3; % range of volume (km3)
range_radius     = 1000.*(range_volume_km3./(4*pi/3)).^(1/3); % range of radius (m)
range_water    = range_water./100; % convert to mass fraction
range_co2      = range_co2./1e6;   % convert to mass fraction
range_vfr        = 10.^range_log_vfr; % volume flow rate (km3/yr)
range_mfr        = param.rho_m0.*range_vfr.*1e9/(3600*24*365); % convert to kg/s

[grid_radius,grid_water,grid_co2,grid_mfr,grid_depth] = ndgrid(range_radius,range_water,range_co2,range_mfr,range_depth);

% number of runs
number_runs      = length(grid_depth(:));
disp(['Number of models to be run = ' num2str(number_runs)])

if save_output == 1
    foldername = [param.composition '_H2O_' num2str(range_water*100) '_CO2_' num2str(range_co2*1e6)];
    mkdir(['output/' foldername])
    mkdir(['figures/' foldername]);
end


% run chamber model over the grid values
for run_s = 1:50
    run_samples=1:18:900; 
    run_i=run_samples(run_s);
    run_n       = run_i+1000;     % Number for input/output files


    % chamber volume and radius 
    a             = grid_radius(run_i);%1000*(10/(4*pi/3))^(1/3);       % initial radius of the chamber (m)
    V_0           = 4*pi/3*a^3; % initial volume of the chamber (m^3)
    
    % volatile concentrations
    InitialConc_H2O = grid_water(run_i); % 450 ppm and 4.2 wt % is interesting as a failure
    InitialConc_CO2 = grid_co2(run_i);
    
    % recharge rate
    mdot_in       = grid_mfr(run_i);     % Mass inflow rate (kg/s)
     
    % depth of reservoir
    depth        = grid_depth(run_i); % depth (m)   
    
    % Run the model
    disp('********************************')
    disp(['Running ' param.composition ' - CO2 version of the code - m_co2 = ' num2str(InitialConc_CO2) ' ; m_h2o = ' num2str(InitialConc_H2O) ' ; depth = ' num2str(depth) ' ; recharge rate = ' num2str(mdot_in)])
    disp(['Running simulation ' num2str(run_i) ' out of ' num2str(number_runs)])
    
    run mainChamber.m
    
end

