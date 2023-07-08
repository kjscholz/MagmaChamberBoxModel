
save_compare=1;
compare_vec=1:18:900;
for compare_count=1:50
    
run_number= 1000+compare_vec(compare_count);
%run_number = compare_vec(compare_count);
composition_folder = 'silicic_H2O_6_CO2_1000/';

% Read in original data
filepath = ['D:/CO2_Paper_Silicic_noDike_highRES/output/' composition_folder 'run_' num2str(run_number) '.mat'];
load(filepath)
original_time = time;
original_P = P;
original_X_co2 = X_co2;
original_eps_g = eps_g;
original_T = T;
original_eps_x = eps_x;
original_tot_Mass = tot_Mass;
original_tot_Mass_0 = tot_Mass(1);
original_tot_Mass_H2O = tot_Mass_H2O;
original_tot_Mass_H2O_0 = tot_Mass_H2O(1);
original_tot_Mass_CO2 = tot_Mass_CO2;
original_tot_Mass_CO2_0 = tot_Mass_CO2(1);
original_V_0=V_0
original_mdot_in = mdot_in

% Read in new data
filepath = ['../usr/output/' composition_folder 'run_' num2str(run_number) '.mat'];
load(filepath)
new_time = time;
new_P = P;
new_X_co2 = X_co2;
new_eps_g = eps_g;
new_T = T;
new_eps_x = eps_x;
new_tot_Mass = tot_Mass;
new_tot_Mass_0 = tot_Mass(1);
new_tot_Mass_H2O = tot_Mass_H2O;
new_tot_Mass_H2O_0 = tot_Mass_H2O(1);
new_tot_Mass_CO2 = tot_Mass_CO2;
new_tot_Mass_CO2_0 = tot_Mass_CO2(1);
% 



% Read in new data
filepath = ['D:/CO2model_newICFinder/MagmaChamberModel_main/usr/output/' composition_folder 'run_' num2str(run_number) '.mat'];
load(filepath)
IC_time = time;
IC_P = P;
IC_X_co2 = X_co2;
IC_eps_g = eps_g;
IC_T = T;
IC_eps_x = eps_x;
IC_tot_Mass = tot_Mass;
IC_tot_Mass_0 = tot_Mass(1);
IC_tot_Mass_H2O = tot_Mass_H2O;
IC_tot_Mass_H2O_0 = tot_Mass_H2O(1);
IC_tot_Mass_CO2 = tot_Mass_CO2;
IC_tot_Mass_CO2_0 = tot_Mass_CO2(1);
% 
% 
p = 6; % point size
ftsize = 12;
fig_size_x = 38;
fig_size_y = 19;

s2yr = 1/(3600*24*365*1e3);
legend_values = ["old"; "new IC, fixed"; "new IC, old enthalpy"];

n=4;

figure('Renderer', 'painters','units','centimeters', 'Position', [10 10 fig_size_x fig_size_y]) 

%t = tiledlayout(n,2,'TileSpacing','Compact','Padding','Compact');

% Plot XCo2 
subplot(2,4,1)

plot(original_time.*s2yr,original_X_co2,'LineWidth',1)
hold on
plot(new_time.*s2yr,new_X_co2,'LineWidth',1)
plot(IC_time.*s2yr,IC_X_co2,'LineWidth',1)
ylabel('X_{co2}')
set(gca,'LineWidth',1, 'FontSize',ftsize)
legend(legend_values)


% eps_g 
subplot(2,4,2)

plot(original_time.*s2yr,original_eps_g,'LineWidth',1)
hold on
plot(new_time.*s2yr,new_eps_g,'LineWidth',1)
plot(IC_time.*s2yr,IC_eps_g,'LineWidth',1)
%xlabel('Time (years)')
ylabel('\epsilon_{g}')
set(gca,'LineWidth',1, 'FontSize',ftsize)





% Plot T 
subplot(2,4,3)

plot(original_time.*s2yr,original_T-273,'LineWidth',1)
hold on
plot(new_time.*s2yr,new_T-273,'LineWidth',1)
plot(IC_time.*s2yr,IC_T-273,'LineWidth',1)
ylabel('T (C)')

set(gca,'LineWidth',1, 'FontSize',ftsize)

% eps_x 
subplot(2,4,4)

plot(original_time.*s2yr,original_eps_x,'LineWidth',1)
hold on
plot(new_time.*s2yr,new_eps_x,'LineWidth',1)
plot(IC_time.*s2yr,IC_eps_x,'LineWidth',1)

%xlabel('Time (years)')
ylabel('\epsilon_{x}')
set(gca,'LineWidth',1, 'FontSize',ftsize)
ylim([0 0.5])


% Plot M 
subplot(2,4,5)
plot(original_time.*s2yr,original_tot_Mass,'LineWidth',1)
hold on
plot(new_time.*s2yr,new_tot_Mass,'LineWidth',1)
plot(IC_time.*s2yr,IC_tot_Mass,'LineWidth',1)
ylabel('M')
xlabel('Time (kyr)')
ylabel('M')
set(gca,'LineWidth',1, 'FontSize',ftsize)


% Plot M_h2o / M_0
subplot(2,4,6)
plot(original_time.*s2yr,original_tot_Mass_H2O./original_tot_Mass_0,'LineWidth',1)
hold on
plot(new_time.*s2yr,new_tot_Mass_H2O./new_tot_Mass_0,'LineWidth',1)
plot(IC_time.*s2yr,IC_tot_Mass_H2O./IC_tot_Mass_0,'LineWidth',1)
ylabel('M_w/M_0')
xlabel('Time (kyr)')
ylabel('M_w/M_0')
set(gca,'LineWidth',1, 'FontSize',ftsize)

% Plot M_co2 / M_0
subplot(2,4,7)
plot(original_time.*s2yr,original_tot_Mass_CO2./original_tot_Mass_0,'LineWidth',1)
hold on
plot(new_time.*s2yr,new_tot_Mass_CO2./new_tot_Mass_0,'LineWidth',1)
plot(IC_time.*s2yr,IC_tot_Mass_CO2./IC_tot_Mass_0,'LineWidth',1)
ylabel('M_c/M_0')
xlabel('Time (kyr)')
ylabel('M_c/M_0')
set(gca,'LineWidth',1, 'FontSize',ftsize)

 %Plot P
subplot(2,4,8)

plot(original_time.*s2yr,original_P/1e6,'LineWidth',1)
hold on
plot(new_time.*s2yr,new_P/1e6,'LineWidth',1)

plot(IC_time.*s2yr,IC_P/1e6,'LineWidth',1)
ylabel('P (MPa)')
xlabel('Time (kyr)')


if save_compare ==1
    mkdir(['comparison_figures/' composition_folder])
    saveas(gca,['comparison_figures/' composition_folder 'run_' num2str(run_number) '.jpg'] )
end
clearvars -except save_compare compare_vec compare_count
close all
end