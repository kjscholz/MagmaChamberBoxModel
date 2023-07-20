save_timeseries=1;
run_list = 1:18:900;
addpath('../src/')
% set up to compare up to three different runs, length of this vector
% determine how many are plotted

%composition_folders = ["silicic_H2O_4_CO2_100/";"silicic_H2O_4_CO2_1000/";"silicic_H2O_6_CO2_100/"];
composition_folders = ["mafic_H2O_0.5_CO2_500/";"mafic_H2O_0.5_CO2_10000/";"mafic_H2O_2_CO2_500/"];

% specify legend entries from composition folders
old = "_";
new = " ";
legend_entries = replace(composition_folders,old,new);

figure_directory = 'figures_timeseries/';

%figure_name_prefix = 'silicic';
figure_name_prefix = 'mafic';

figure_type= '.jpg'; %options include .jpg and .epsc

for count = 1:3%length(run_list)   
run_prefix = 1000;
run_number = run_prefix + run_list(count);

% Read in first run
filepath = ['../usr/output/' convertStringsToChars(composition_folders(1)) 'run_' num2str(run_number) '.mat'];
load(filepath)
time_1 = time;
P_1 = P;
X_co2_1 = X_co2;
eps_g_1 = eps_g;
T_1 = T;
eps_x_1 = eps_x;
tot_Mass_1 = tot_Mass;
tot_Mass_ini_1 = tot_Mass(1);
tot_Mass_H2O_1 = tot_Mass_H2O;
tot_Mass_H2O_ini_1 = tot_Mass_H2O(1);
m_H2O_diss_1 = M_h2o_diss;
tot_Mass_CO2_1 = tot_Mass_CO2;
m_CO2_diss_1 = M_co2_diss;
tot_Mass_CO2_ini_1 = tot_Mass_CO2(1);
V0 = V_0
mdot_in = mdot_in

if length(composition_folders)>1
% Read in second run
filepath = ['../usr/output/' convertStringsToChars(composition_folders(2)) 'run_' num2str(run_number) '.mat'];
load(filepath)
time_2 = time;
P_2 = P;
X_co2_2 = X_co2;
eps_g_2 = eps_g;
T_2 = T;
eps_x_2 = eps_x;
tot_Mass_2 = tot_Mass;
tot_Mass_ini_2 = tot_Mass(1);
tot_Mass_H2O_2 = tot_Mass_H2O;
m_H2O_diss_2 = M_h2o_diss;
tot_Mass_H2O_ini_2 = tot_Mass_H2O(1);
tot_Mass_CO2_2 = tot_Mass_CO2;
m_CO2_diss_2 = M_co2_diss;
tot_Mass_CO2_ini_2 = tot_Mass_CO2(1);
% 

end

if length(composition_folders)>2
% Read third run
filepath = ['../usr/output/' convertStringsToChars(composition_folders(3)) 'run_' num2str(run_number) '.mat'];
load(filepath)
time_3 = time;
P_3 = P;
X_co2_3 = X_co2;
eps_g_3 = eps_g;
T_3 = T;
eps_x_3 = eps_x;
tot_Mass_3 = tot_Mass;
tot_Mass_ini_3 = tot_Mass(1);
tot_Mass_H2O_3 = tot_Mass_H2O;
tot_Mass_H2O_ini_3 = tot_Mass_H2O(1);
m_H2O_diss_3 = M_h2o_diss;
tot_Mass_CO2_3 = tot_Mass_CO2;
m_CO2_diss_3 = M_co2_diss;
tot_Mass_CO2_ini_3 = tot_Mass_CO2(1);
end
% 
% 
p = 6; % point size
ftsize = 12;
fig_size_x = 38;
fig_size_y = 19;

s2yr = 1/(3600*24*365*1e3);



figure('Renderer', 'painters','units','centimeters', 'Position', [10 10 fig_size_x fig_size_y]) 

% Plot XCo2 
subplot(2,4,1)

plot(time_1.*s2yr,X_co2_1,'LineWidth',1)
hold on

if length(composition_folders)>1
plot(time_2.*s2yr,X_co2_2,'LineWidth',1)
end

if length(composition_folders)>2
plot(time_3.*s2yr,X_co2_3,'LineWidth',1)
end

ylabel('X_{co2}')
set(gca,'LineWidth',1, 'FontSize',ftsize)
legend(legend_entries)


% eps_g 
subplot(2,4,2)

plot(time_1.*s2yr,eps_g_1,'LineWidth',1)
hold on

if length(composition_folders)>1
plot(time_2.*s2yr,eps_g_2,'LineWidth',1)
end

if length(composition_folders)>2
plot(time_3.*s2yr,eps_g_3,'LineWidth',1)
end

%xlabel('Time (years)')
ylabel('\epsilon_{g}')
set(gca,'LineWidth',1, 'FontSize',ftsize)


% Plot T 
subplot(2,4,3)

plot(time_1.*s2yr,T_1-273,'LineWidth',1)
hold on
if length(composition_folders)>1
plot(time_2.*s2yr,T_2-273,'LineWidth',1)
end
if length(composition_folders)>2
plot(time_3.*s2yr,T_3-273,'LineWidth',1)
end
ylabel('T (C)')

set(gca,'LineWidth',1, 'FontSize',ftsize)

% eps_x 
subplot(2,4,4)

plot(time_1.*s2yr,eps_x_1,'LineWidth',1)
hold on

if length(composition_folders)>1
plot(time_2.*s2yr,eps_x_2,'LineWidth',1)
end

if length(composition_folders)>2
plot(time_3.*s2yr,eps_x_3,'LineWidth',1)
end

%xlabel('Time (years)')
ylabel('\epsilon_{x}')
set(gca,'LineWidth',1, 'FontSize',ftsize)
ylim([0 0.5])


% Plot M 
subplot(2,4,5)
plot(time_1.*s2yr,tot_Mass_1./tot_Mass_ini_1,'LineWidth',1)
hold on

if length(composition_folders)>1
plot(time_2.*s2yr,tot_Mass_2./tot_Mass_ini_2,'LineWidth',1)
end

if length(composition_folders)>2
plot(time_3.*s2yr,tot_Mass_3./tot_Mass_ini_3,'LineWidth',1)
end 

ylabel('M/M_0')
xlabel('Time (kyr)')

set(gca,'LineWidth',1, 'FontSize',ftsize)


% Plot M_h2o / M_0
subplot(2,4,6)
plot(time_1.*s2yr,tot_Mass_H2O_1./tot_Mass_1,'Color',[0 0.4470 0.7410],'LineStyle','- -','LineWidth',1)
hold on
plot(time_1.*s2yr,m_H2O_diss_1./100,'Color',[0 0.4470 0.7410],'LineWidth',1)

if length(composition_folders)>1
plot(time_2.*s2yr,tot_Mass_H2O_2./tot_Mass_2,'Color',[0.8500 0.3250 0.0980],'LineStyle','- -','LineWidth',1)
plot(time_2.*s2yr,m_H2O_diss_2./100,'Color',[0.8500 0.3250 0.0980],'LineWidth',1)
end

if length(composition_folders)>2
plot(time_3.*s2yr,tot_Mass_H2O_3./tot_Mass_3,'Color',[0.9290 0.6940 0.1250],'LineStyle','- -','LineWidth',1)
plot(time_3.*s2yr,m_H2O_diss_3./100,'Color',[0.9290 0.6940 0.1250],'LineWidth',1)
end

ylabel('H_2O Concentration')
xlabel('Time (kyr)')
legend("total", "in melt")
set(gca,'LineWidth',1, 'FontSize',ftsize)

% Plot M_co2 / M_0
subplot(2,4,7)
plot(time_1.*s2yr,tot_Mass_CO2_1./tot_Mass_1,'Color',[0 0.4470 0.7410],'LineStyle','- -','LineWidth',1)
hold on
plot(time_1.*s2yr,m_CO2_diss_1./1e6,'Color',[0 0.4470 0.7410],'LineWidth',1)

if length(composition_folders)>1
plot(time_2.*s2yr,tot_Mass_CO2_2./tot_Mass_2,'Color',[0.8500 0.3250 0.0980],'LineStyle','- -','LineWidth',1)
plot(time_2.*s2yr,m_CO2_diss_2./1e6,'Color',[0.8500 0.3250 0.0980],'LineWidth',1)
end

if length(composition_folders)>2
plot(time_3.*s2yr,tot_Mass_CO2_3./tot_Mass_3,'Color',[0.9290 0.6940 0.1250],'LineStyle','- -','LineWidth',1)
plot(time_3.*s2yr,m_CO2_diss_3./1e6,'Color',[0.9290 0.6940 0.1250],'LineWidth',1)
end
legend("total", "in melt")
ylabel('CO_2 Concentration')
xlabel('Time (kyr)')

set(gca,'LineWidth',1, 'FontSize',ftsize)

% Plot P
subplot(2,4,8)

plot(time_1.*s2yr,P_1/1e6,'LineWidth',1)
hold on

if length(composition_folders)>1
plot(time_2.*s2yr,P_2/1e6,'LineWidth',1)
end

if length(composition_folders)>2
plot(time_3.*s2yr,P_3/1e6,'LineWidth',1)
end
ylabel('P (MPa)')
xlabel('Time (kyr)')

sgtitle(['Run Number ', num2str(run_list(count)), ', V_0 = ', num2str(round((V0/1e9),3,'significant')), ' km^3, Mdot_{in} = ',num2str(round((mdot_in/1e9),4,'significant')), ' kg/s']);

if save_timeseries ==1
    mkdir(figure_directory)
    saveas(gca,[figure_directory, figure_name_prefix, '_run_', num2str(run_number), figure_type] )
    close 
end
clearvars -except save_timeseries figure_directory figure_name_prefix figure_type run_list composition_folders legend_entries

end