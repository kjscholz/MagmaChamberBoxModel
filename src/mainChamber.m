tic

global storeTime storeTemp Mdot_in_pass Mdot_out_pass P_lit
global phase tot_h2o_frac_in tot_co2_frac_in
global switch_Tprofile Q_out_old lengthTime maxTime

% initial pressure
P_lit    = param.rho_r*grav_acc*depth; % lithostatic pressure (Pa)
P_0      = P_lit; % initial chamber pressure (Pa)

% initial chamber temperature (based on initial crystal fraction)
if strcmp(param.composition,'silicic')
    T_0   = find_liq_silicic(InitialConc_H2O,InitialConc_CO2,P_0,ini_eps_x);
elseif strcmp(param.composition,'mafic')
    T_0   = find_liq_mafic(InitialConc_H2O,InitialConc_CO2,P_0,ini_eps_x);
end

% set the mass inflow rate and conditions
Mdot_in_pass  = mdot_in;
Mdot_out_pass = 10000;
if strcmp(T_R,'same')
    T_R           = T_0;       % Temperature of recharging magma (K)
end
T_in = T_R; % set temperature of recharging magma

% far-field temperature
param.Tb = T_surface+T_gradient*depth; % background temperature crust (K)

lengthTime     = 0;
maxTime        = 0;
switch_Tprofile= 0;
Q_out_old      = 0;

% time
begin_time     = 0;         % initialize time
end_time       = max_simulation_time;      % maximum simulation time in seconds

rho_m0 = param.rho_m0;      % initial melt density (kg/m^3)
rho_x0 = param.rho_x0;      % initial melt density (kg/m^3)

[rho_g0,~,~]  = eos_g(P_0,T_0); % initial gas density

if strcmp(param.composition,'silicic')
    [eps_x0,~,~,~,~,~]=crystal_fraction_silicic(T_0,P_0,0,0,InitialConc_H2O,InitialConc_CO2,rho_m0,rho_x0,rho_g0);
elseif strcmp(param.composition,'mafic')
    [eps_x0,~,~,~,~,~]=crystal_fraction_mafic(T_0,P_0,0,0,InitialConc_H2O,InitialConc_CO2,rho_m0,rho_x0,rho_g0);
end

eps_m0=1-eps_x0;
rho = rho_m0*eps_m0+rho_x0*eps_x0;
M_tot =  V_0*rho; % Total mass, initial
M_co2_0 = V_0*rho*InitialConc_CO2; % Total mass of CO2, initial
M_h2o_0 = InitialConc_H2O*V_0*rho; % Total mass of CO2, initial

%if InitialConc_H2O > m_eq0
if strcmp(param.composition,'silicic')
    [eps_g0, X_co20, ~]=IC_Finder_silicic(M_h2o_0, M_co2_0, M_tot, P_0, T_0, V_0, rho_m0,rho_x0);
elseif strcmp(param.composition,'mafic')
    [eps_g0, X_co20, ~]=IC_Finder_mafic(M_h2o_0, M_co2_0, M_tot, P_0, T_0, V_0, rho_m0,rho_x0);
end
disp(['1st ICfinder done, starting in ' num2str(phase) '-phase, eps_g0 = ' num2str(eps_g0)])


% update initial crystal volume fraction
if strcmp(param.composition,'silicic')
    [eps_x0,~,~,~]=crystal_fraction_silicic(T_0,P_0,eps_g0,X_co20,InitialConc_H2O,InitialConc_CO2,rho_m0,rho_x0,rho_g0);%crystal_fraction(T,P, eps_g,X_co2,mH20,mCO2,rho_m,rho_x,rho_g)
elseif strcmp(param.composition,'mafic')
    [eps_x0,~,~,~]=crystal_fraction_mafic(T_0,P_0,eps_g0,X_co20,InitialConc_H2O,InitialConc_CO2,rho_m0,rho_x0,rho_g0);%crystal_fraction(T,P, eps_g,X_co2,mH20,mCO2,rho_m,rho_x,rho_g)
end

% update initial melt volume fraction
eps_m0 = 1-eps_x0-eps_g0;
% update initial bulk density (kg/m^3)
rho_0  = (1-eps_g0-eps_x0)*rho_m0 + eps_g0*rho_g0 + eps_x0*rho_x0;

% update solubility
if phase == 2
    m_eq0 = M_h2o_0/((1-eps_x0)*rho_m0*V_0);
    C_co20 = M_co2_0/((1-eps_x0)*rho_m0*V_0);
    X_co20 = 0;
    test_phase=1;
else
    if strcmp(param.composition,'silicic')
        [m_eq0,~,~,~,C_co20,~,~,~]   = exsolve_silicic(P_0,T_0,X_co20);
    elseif strcmp(param.composition,'mafic')
        [m_eq0,~,~,~,C_co20,~,~,~]   = exsolve_mafic(P_0,T_0,X_co20);
    end
end

% Calculate the total water content for initial magma
m_h2o_tot  = M_h2o_0;

% Calculate the total CO2 content for initial magma
m_co2_tot  = M_co2_0;

% Calculate the water content (concentration) for inflowing magma contents
tot_h2o_frac_in = M_h2o_0/(rho_0*V_0);
tot_co2_frac_in = M_co2_0/(rho_0*V_0);

tot_Mass_0 = V_0*rho_0;
tot_Mass_H2O_0=M_h2o_0;
tot_Mass_CO2_0=M_co2_0;


sw.heat_cond     = 1; % switch cooling module on/off
sw.visc_relax    = 1; % switch viscous relaxation on/off

% initialize vector to store quantities
storeTime = [];
storeTemp = [];
time  = [];
P     = [];
T     = [];
eps_g = [];
eps_g_test=[];
X_co2_test=[];
eps_x = [];
deps_x_dP = [];
V     = [];
rho_m = [];
rho_x = [];
rho_g = [];
tot_Mass = [];
tot_Mass_H2O = [];
tot_Mass_CO2 = [];
X_co2 = [];
time_erupt = [];
dur_erupt = [];
mass_erupt = [];
vol_erupt = [];
epsg_erupt = [];
epsx_erupt = [];
mass_added = [];
dur_cycle = [];
er_end_times = [];

ii      = 1; % initial iteration for while loop
ii_max  = max_num_eruptions; % maximum iteration in while loop
phase0  = 0; % number that indicates why the code stopped

% Initial temperature and viscosity profile around chamber
[quadpts,weights] = GLQ_points_weights_hard(param.GLQ_n);
cc    = 10*a; %outer radius for heat conduction (m)
b     = a +cc;
quadpts_r = (b-a)/2.*quadpts + (a+b)/2;
storeTime =0;
storeTemp = T_0;
Trt = heat_conduction_chamber_profile(param.maxn,a,cc,quadpts_r,param.kappa,param.Tb);
A = param.A; % material-dependent constant for viscosity law (Pa s)
B = param.B; % molar gas constant (J/mol/K)
G = param.G; % activation energy for creep (J/mol)
eta_rt     = A.*exp(G./B./Trt);
I          = (b-a)/2.*sum(weights.*(eta_rt./(quadpts_r).^4));
eta_r0      = 3*a^3*I;
disp('crust viscosity (Pa s)')
disp(eta_r0)

% time scales (s)
tau_inj            = rho_0*V_0/(mdot_in); % injection timescale
tau_cooling        = a^2/param.kappa;     % cooling timescale
tau_visco_elastic  = eta_r0/P_crit;       % relaxation timescale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



sw.eruption = 0;

while begin_time < end_time && phase0 ==0  && ii<ii_max
    % set up ode
    tspan    = [begin_time  end_time];  % solve from t=0s to t=end_time
    IC       = [P_0 T_0 eps_g0 V_0 rho_m0 rho_x0 X_co20 tot_Mass_0 tot_Mass_H2O_0 tot_Mass_CO2_0]; % store initial conditions
    
    options     = odeset('RelTol',rel_tol,'AbsTol',abs_tol,...  % used to be 1e7, 1
        'Events',@(z,y) stopChamber(y,sw.eruption, P_lit,P_crit,param));%'Refine',1,'MaxStep',1e7,'InitialStep',1e3);
    % solve ode
    [X,Y,TE,YE,IE] = ode15s(@(t,y) odeChamber(t,y,        ...
        sw, ...
        T_in,...
        param), ...
        tspan,    ...
        IC,       ...
        options);
    
    % store output, not the most proper way of coding but it works
    time  = [time;  X];
    P     = [P;     Y(:,1)];
    T     = [T;     Y(:,2)];
    eps_g = [eps_g; Y(:,3)];
    V     = [V;     Y(:,4)];
    rho_m = [rho_m; Y(:,5)];
    rho_x = [rho_x; Y(:,6)];
    X_co2 = [X_co2; Y(:,7)];
    tot_Mass = [tot_Mass; Y(:,8)];
    tot_Mass_H2O = [tot_Mass_H2O; Y(:,9)];
    tot_Mass_CO2 = [tot_Mass_CO2; Y(:,10)];
    
    % re-initialize
    begin_time = time(end);
    P_0        = P(end);
    T_0        = T(end);
    eps_g0     = eps_g(end);
    V_0        = V(end);
    rho_m0     = rho_m(end);
    rho_x0     = rho_x(end);
    X_co20     = X_co2(end);
    tot_Mass_0 = tot_Mass(end);
    tot_Mass_H2O_0 = tot_Mass_H2O(end);
    tot_Mass_CO2_0 = tot_Mass_CO2(end);
    M_h2o_0 =tot_Mass_H2O_0;
    M_co2_0 = tot_Mass_CO2_0;
    M_tot = tot_Mass_0;
    storeTime = storeTime(storeTime<time(end));
    storeTemp = storeTemp(storeTime<time(end));
    
    m_h2o=tot_Mass_H2O(end)/tot_Mass(end);
    m_co2=tot_Mass_CO2(end)/tot_Mass(end);
    [rho_g0,~,~]  = eos_g(P_0,T_0); % initial gas density
    
    if strcmp(param.composition,'silicic')
        [eps_x0,~,~,~] =  crystal_fraction_silicic(T_0,P_0,eps_g0,X_co20,m_h2o,m_co2,rho_m0,rho_x0,rho_g0);
    elseif strcmp(param.composition,'mafic')
        [eps_x0,~,~,~] =  crystal_fraction_mafic(T_0,P_0,eps_g0,X_co20,m_h2o,m_co2,rho_m0,rho_x0,rho_g0);
    end
    
    % If it just reached critical pressure and need to start an eruption
    if ismember(3,IE)  && eps_x0<0.5 % && ~ismember(7,IE)
        er_start_time = time(end);
        time_erupt = [time_erupt er_start_time];
        sw.eruption=1;
        er_t0 = time(end);
        rho_er_start = (1-eps_g0-eps_x0).*rho_m0 + eps_g0.*rho_g0 + eps_x0.*rho_x0;
        
        epsg_erupt  = [epsg_erupt eps_g0];
        epsx_erupt  = [epsx_erupt eps_x0];
        
    end
    
    % If it just finished an eruption...
    if ismember(4,IE)
        % store eruption times and durations
        er_end_time = time(end);
        dur = er_end_time-er_start_time;
        mass_withdrawn = dur*Mdot_out_pass;
        vol_withdrawn = mass_withdrawn/rho_er_start;
        
        er_ind = find(time==er_start_time)+1;
        
        % store eruption times and durations
        dur_erupt   = [dur_erupt dur];
        mass_erupt  = [mass_erupt mass_withdrawn];
        vol_erupt   = [vol_erupt vol_withdrawn];
        er_end_times = [er_end_times er_end_time];
        
        k = length(time_erupt);
        if k > 1
            dur_thiscycle = time(end)-er_end_times(k-1);
        else
            dur_thiscycle = time(end);
        end
        dur_cycle = [dur_cycle dur_thiscycle];
        mass_added = [mass_added dur_thiscycle*mdot_in];
        
        sw.eruption = 0;
        ii = ii+1; % add to the eruption count
    end
    
    %If it needs to transition between 2 and 3 phases
    if  ismember(6,IE) || ismember(8,IE)
        
        phase_here=phase;
        
        if strcmp(param.composition,'silicic')
            [eps_g_temp, X_co2_temp, C_co2_temp]=IC_Finder_silicic(tot_Mass_H2O_0,...
                tot_Mass_CO2_0, tot_Mass_0, P_0, T_0, V_0, rho_m0,rho_x0);
        elseif strcmp(param.composition,'mafic')
            [eps_g_temp, X_co2_temp, C_co2_temp]=IC_Finder_mafic(tot_Mass_H2O_0,...
                tot_Mass_CO2_0, tot_Mass_0, P_0, T_0, V_0, rho_m0,rho_x0);
        end
        
        if phase_here ~= phase
            eps_g0 = eps_g_temp;
            X_co20 = X_co2_temp;
            C_co2  = C_co2_temp;
            disp(['went from ' num2str(phase_here) '-phase to ' num2str(phase) '-phase'])
        end
        
        
    end
    
    
    if IE==1
        disp('eps_x became 0.')
        phase0=1;
    elseif IE==2
        disp('eps_x/(1-eps_g) became 0.8')
        phase0=2;
    elseif ismember(5,IE)
        disp('eps_x became 0.5')
        phase0=5;
    elseif ismember(7,IE)
        disp('too much underpressure - collapse')
        phase0=7;
    elseif isempty(IE)
        disp('you reached the end of time')
        phase0=8;
    end
    
end





%% figures and post process calculations
%
% %crystal volume fraction
eps_x = zeros(size(time));
rho_g = zeros(size(time));
deps_x_dP = zeros(size(time));
mh2o=tot_Mass_H2O./tot_Mass;
mco2=tot_Mass_CO2./tot_Mass;
for i = 1:length(time)
    [rho_g(i),~,~]  = eos_g(P(i),T(i));
    
    if strcmp(param.composition,'silicic')
        [eps_x(i),deps_x_dP(i),~,~] =  crystal_fraction_silicic(T(i),P(i),eps_g(i),X_co2(i),mh2o(i),mco2(i),rho_m(i),rho_x(i),rho_g(i));
        
    elseif strcmp(param.composition,'mafic')
        [eps_x(i),deps_x_dP(i),~,~] =  crystal_fraction_mafic(T(i),P(i),eps_g(i),X_co2(i),mh2o(i),mco2(i),rho_m(i),rho_x(i),rho_g(i));
        
    end
    
end

%
% % dissolved water mass fraction
m_eq = zeros(size(P));
dmeqdP = zeros(size(P));
dmeqdXco2= zeros(size(P));
dC_co2dP= zeros(size(P));
dC_co2dXco2= zeros(size(P));
C_co2    = zeros(size(P));

for i = 1:length(m_eq)
    if strcmp(param.composition,'silicic')
        [m_eq(i),dmeqdP(i),~,dmeqdXco2(i),C_co2(i),dC_co2dP(i),~,dC_co2dXco2(i)] = exsolve_silicic(P(i),T(i), X_co2(i));
        
    elseif strcmp(param.composition,'mafic')
        [m_eq(i),dmeqdP(i),~,dmeqdXco2(i),C_co2(i),dC_co2dP(i),~,dC_co2dXco2(i)] = exsolve_mafic(P(i),T(i), X_co2(i));
        
    end
end
%
% % gas density
rho_g = zeros(size(P));
drho_g_dP = zeros(size(P));
for i = 1:length(rho_g)
    [rho_g(i),drho_g_dP(i),~] = eos_g(P(i),T(i));
end
%
% % bulk density
rho  = (1-eps_g-eps_x).*rho_m + eps_g.*rho_g + eps_x.*rho_x;

% % bulk heat capacity
[c_g, ~]=gas_heat_capacity(X_co2);
c  = ((1-eps_g-eps_x).*rho_m.*param.c_m + eps_g.*rho_g.*c_g + eps_x.*rho_x.*param.c_x)./rho;



% compute mass of co2 and h2o dissolved
 % 1. first get mass of co2 and h2o in gas
    M_co2_gas = X_co2.*param.mm_co2.*(rho_g.*eps_g.*V./(X_co2.*param.mm_co2+...
        (1-X_co2).*param.mm_h2o));
    M_h2o_gas = (1-X_co2).*param.mm_h2o.*(rho_g.*eps_g.*V./(X_co2.*param.mm_co2+...
        (1-X_co2).*param.mm_h2o));
    % 2. subtract from total mass and convert to wt% or ppm
    M_co2_diss = 1e6*(tot_Mass_CO2-M_co2_gas)./(rho_m.*V.*(1-eps_g-eps_x)); % ppm
    M_h2o_diss = 100*(tot_Mass_H2O-M_h2o_gas)./(rho_m.*V.*(1-eps_g-eps_x)); % wt%
    




time_of_run = toc

if save_output == 1
    save(['output/' foldername '/run_', num2str(run_n), '.mat'])
    
end

if save_figures == 1
    figure('Renderer', 'painters', 'Position', [10 10 900 900])
    
    subplot(5,2,1)
    plot(time./(3600*24*365),T-273,'LineWidth',2)
    xlabel('time (yr)','FontSize',16), ylabel('T','FontSize',16)
    set(gca,'FontSize',12)
    
    subplot(5,2,2)
    plot(time./(3600*24*365),eps_x,'LineWidth',2)
    xlabel('time (yr)','FontSize',16), ylabel('eps_x','FontSize',16)
    set(gca,'FontSize',12)
    
    subplot(5,2,3)
    plot(time./(3600*24*365),eps_g,'LineWidth',2)
    xlabel('time (yr)','FontSize',16), ylabel('eps_g','FontSize',16)
    set(gca,'FontSize',12)
    
    subplot(5,2,4)
    plot(time./(3600*24*365),P./1e6,'LineWidth',2)
    xlabel('time (yr)','FontSize',16), ylabel('P(MPa)','FontSize',16)
    set(gca,'FontSize',12)
    
    subplot(5,2,5)
    plot(time./(3600*24*365),X_co2,'LineWidth',2)
    xlabel('time (yr)','FontSize',16), ylabel('X_co2','FontSize',16)
    set(gca,'FontSize',12)
    
    subplot(5,2,6)
    plot(time./(3600*24*365),m_eq,'LineWidth',2)
    xlabel('time (yr)','FontSize',16), ylabel('m_eq','FontSize',16)
    set(gca,'FontSize',12)
    
    subplot(5,2,7)
    plot(time./(3600*24*365),C_co2,'LineWidth',2)
    xlabel('time (yr)','FontSize',16), ylabel( 'C_co2','FontSize',16)
    set(gca,'FontSize',12)
    
    subplot(5,2,8)
    plot(time./(3600*24*365),tot_Mass,'LineWidth',2)
    xlabel('time (yr)','FontSize',16), ylabel('total Mass solved','FontSize',16)
    set(gca,'FontSize',12)
    
    subplot(5,2,9)
    plot(time./(3600*24*365),rho.*V,'LineWidth',2)
    xlabel('time (yr)','FontSize',16), ylabel('rho * V (total mass)','FontSize',16)
    set(gca,'FontSize',12)
    
    subplot(5,2,10)
    plot(time./(3600*24*365),V./1e9,'LineWidth',2)
    xlabel('time (yr)','FontSize',16), ylabel('Volume (km^3)','FontSize',16)
    set(gca,'FontSize',12)
    
    sgtitle([param.composition ', H2O ' num2str(range_water*100) ', CO2 ' num2str(range_co2*1e6) ', run ' num2str(run_n) ', Mdot_{in} ' num2str(mdot_in) ' kg/s, V ' num2str(V(1)/1e9) ' km^3' ])
    figname=['figures/' foldername '/run_' num2str(run_n) '.jpg'];
    
    saveas(gcf, figname)
    close all
end

