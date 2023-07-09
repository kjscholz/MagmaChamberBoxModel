% Use this to run one model at a time and display output
clear all
close all
clc
addpath('../src/')

run_n        = 1; % number for naming input/output files
save_output  = 0; % set to 1 to save output as .mat file
save_figures = 0; % set to 1 to save figures

%% Set initial & boundary conditions

param.composition = 'silicic';
%param.composition = 'mafic';

% chamber volume (km^3) -- note: enter the LOG of volume
V_ch_km3   = 10; % initial volume of chamber (km^3)

% volatile content -- total mass fraction (not just dissolved)
water      = 6;   % wt%
co2        = 1000; % ppm

% recharge rate (km3/yr) -- note: enter the LOG of recharge rate
log_vfr    = -5; % log volume flow rate (km3/yr)

% storage depth (m)
depth      = 8e3;

%% set other properties, do some conversions, and run the code
run set_properties
run set_numerical_params

% conversions
V_0 = V_ch_km3*1e9; % chamber volume in m^3
a   = (V_0*(3/4)*(1/pi))^(1/3); % chamber radius (m)
InitialConc_H2O = water/100; % convert to mass fraction
InitialConc_CO2   = co2/1e6;   % convert to mass fraction
vfr           = 10.^log_vfr; % volume flow rate (km3/yr)
mdot_in        = param.rho_m0.*vfr.*1e9/(3600*24*365); % convert to mass recharge rate (kg/s)


foldername = [param.composition '_H2O_' num2str(water) '_CO2_' num2str(co2)];
if save_output
    mkdir(['output/' foldername]);
end
if save_figures
    mkdir(['figures/' foldername]);
end


% Run the model
disp('********************************')
disp(['Running ' param.composition ' - CO2 version of the code - m_co2 = ' num2str(InitialConc_CO2) ' ; m_h2o = ' num2str(InitialConc_H2O) ' ; depth = ' num2str(depth) ' ; recharge rate = ' num2str(mdot_in)])

run mainChamber.m


%% Plot output
    figure('Renderer', 'painters', 'Position', [0 0 900 900])
    subplot(5,2,1)
    plot(time./(3600*24*365),T-273,'LineWidth',2)
    xlabel('time (yr)','FontSize',16), ylabel('T (C)','FontSize',16)
    set(gca,'FontSize',12)
    
    subplot(5,2,2)
    plot(time./(3600*24*365),eps_x,'LineWidth',2)
    xlabel('time (yr)','FontSize',16), ylabel('\epsilon_x','FontSize',16)
    set(gca,'FontSize',12)
    
    subplot(5,2,3)
    plot(time./(3600*24*365),eps_g,'LineWidth',2)
    xlabel('time (yr)','FontSize',16), ylabel('\epsilon_g','FontSize',16)
    set(gca,'FontSize',12)
    
    subplot(5,2,4)
    plot(time./(3600*24*365),P./1e6,'LineWidth',2)
    xlabel('time (yr)','FontSize',16), ylabel('P (MPa)','FontSize',16)
    set(gca,'FontSize',12)
    
    subplot(5,2,5)
    plot(time./(3600*24*365),X_co2,'LineWidth',2)
    xlabel('time (yr)','FontSize',16), ylabel('X_{co_2}','FontSize',16)
    set(gca,'FontSize',12)
    
    subplot(5,2,6)
    plot(time./(3600*24*365),M_h2o_diss,'LineWidth',2)
    xlabel('time (yr)','FontSize',16), ylabel('dissolved h2o (wt%)','FontSize',16)
    set(gca,'FontSize',12)
    
    subplot(5,2,7)
    plot(time./(3600*24*365),M_co2_diss,'LineWidth',2)
    xlabel('time (yr)','FontSize',16), ylabel( 'dissolved co2 (ppm)','FontSize',16)
    set(gca,'FontSize',12)
    
    subplot(5,2,8)
    plot(time./(3600*24*365),tot_Mass,'LineWidth',2)
    xlabel('time (yr)','FontSize',16), ylabel('total mass (kg)','FontSize',16)
    set(gca,'FontSize',12)
    
    subplot(5,2,9)
    plot(time./(3600*24*365),rho,'LineWidth',2)
    xlabel('time (yr)','FontSize',16), ylabel('magma density (kg/m^3)','FontSize',16)
    set(gca,'FontSize',12)
    
    subplot(5,2,10)
    plot(time./(3600*24*365),V./1e9,'LineWidth',2)
    xlabel('time (yr)','FontSize',16), ylabel('Volume (km^3)','FontSize',16)
    set(gca,'FontSize',12)
    
    sgtitle([param.composition ', H2O = ' num2str(water) ' wt%, CO2 = ' num2str(co2) ' ppm, recharge = ' num2str(10^log_vfr) ' km^3/yr, V = ' num2str(V(1)/1e9) ' km^3' ])
   