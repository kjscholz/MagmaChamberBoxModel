function [eps_g0, X_co20, mco2_diss, eps_x0]=IC_Finder_silicic(M_h2o_0,...
    M_co2_0, M_tot, P_0, T_0,V_0, rho_m0,rho_x0)

% % % 
% NOTE: This version was written by Meredith Townsend on June 6, 2023
% The algorithm is a pseudo MCMC approach to find eps_g0 and X_co20
% % %

global phase
global ss err_tol maxcount_MT
global random_guess eps_g_exp_ini_pass X_co2_ini_pass

param.mm_co2= 44.01e-3; % molar mass of CO2
param.mm_h2o= 18.02e-3; % molar mass of H2O

[rho_g0,~,~]=eos_g(P_0,T_0); % Initial gas density
mH2O=M_h2o_0/M_tot; % total mass fraction of H2O
mCO2=M_co2_0/M_tot; % total mass fraction of CO2
[eps_x0, ~,~, ~]= crystal_fraction_silicic(T_0,P_0, 0,0,mH2O,mCO2,...
    rho_m0,rho_x0,rho_g0); % Initial crystal fraction assume zero gas
eps_m0=1-eps_x0; % Initial melt fraction assuming zero gas

Conc_Water = M_h2o_0/(V_0*rho_m0*eps_m0); % concentration of H2O dissolved in melt if no gas
Conc_co2 = M_co2_0/(V_0*rho_m0*eps_m0); % concentration of CO2 dissolved in melt if no gas

%CHECK IF SATURATED
m_eq_max = exsolve_silicic(P_0, T_0, 0); % Commpute H2O solubility
if Conc_Water >  m_eq_max % If current concentration of H2O exceeds solubility
    phase=3;
else
    [C_co2_sat,~]= exsolve3_silicic(P_0,T_0, Conc_Water); % Commpute CO2 solubility
    if Conc_co2 >  C_co2_sat % If current concentration of CO2 exceeds solubility
        phase = 3;
    else
        phase =2; % If not saturated, set to two-phase
    end
end

if phase ==2
    X_co20=0;
    eps_g0=0;
    mco2_diss=Conc_co2;
    
else

%IF SATURATED, FIND EPS_G0 & X_CO20
    
    % bounds for eps_g and X_co2 (note bounds are on the exponent of eps_g)
    b.eps_g_exp_high = -0.25; b.eps_g_exp_low = -6;
    b.X_co2_high = 0.99999; b.X_co2_low = 0.00001;
    
    % If starting from a random guess, generate that. Otherwise, use what
    % was passed via global variables
    if random_guess==1
        eps_g_exp_ini=b.eps_g_exp_low+(b.eps_g_exp_high-b.eps_g_exp_low)*rand(1,1);
        eps_g_ini = 10^eps_g_exp_ini;
        X_co2_ini = rand(1,1)*(b.X_co2_high-b.X_co2_low)+b.X_co2_low;
    else
        eps_g_exp_ini=eps_g_exp_ini_pass;
        eps_g_ini = 10^eps_g_exp_ini;
        X_co2_ini=X_co2_ini_pass;
    end
    
    % compute predicted mass of co2 and h2p
    [eps_x0, ~,~, ~]= crystal_fraction_silicic(T_0,P_0, eps_g_ini,X_co2_ini,mH2O,mCO2,...
        rho_m0,rho_x0,rho_g0);
    
    eps_m0 = 1-eps_g_ini-eps_x0;
    mm_gas = X_co2_ini*param.mm_co2+(1-X_co2_ini)*param.mm_h2o;
    
    M_co2_predicted = mco2_dissolved_sat(X_co2_ini,P_0,T_0)*eps_m0*V_0*rho_m0+...
        eps_g_ini*rho_g0*V_0*X_co2_ini*param.mm_co2/mm_gas;
    
    M_h2o_predicted = meq_water(X_co2_ini,P_0,T_0)*eps_m0*V_0*rho_m0+...
        eps_g_ini*rho_g0*V_0*(1-X_co2_ini)*param.mm_h2o/mm_gas;
    
    % compute error from initial guess
    err_co2_ini = abs(M_co2_predicted-M_co2_0)/M_co2_0;
    err_h2o_ini = abs(M_h2o_predicted-M_h2o_0)/M_h2o_0;
    err_tot_ini= err_co2_ini+ err_h2o_ini;
    
    count = 0;
    
    while err_tot_ini > err_tol && count < maxcount_MT
        
        count = count+1;
        % take a step to produce a candidate
        eps_g_exp_step = (b.eps_g_exp_high-b.eps_g_exp_low)*ss*randn(1,1);
        eps_g_exp_cand = eps_g_exp_ini+eps_g_exp_step;
        eps_g_cand = 10^eps_g_exp_cand;
        X_co2_cand = X_co2_ini+(b.X_co2_high-b.X_co2_low)*ss*randn(1,1);
        
        % make sure it's in bounds
        while inbounds(eps_g_exp_cand,X_co2_cand,b)
            eps_g_exp_step = (b.eps_g_exp_high-b.eps_g_exp_low)*ss*randn(1,1);
            eps_g_exp_cand = eps_g_exp_ini+eps_g_exp_step;
            eps_g_cand = 10^eps_g_exp_cand;
            X_co2_cand = X_co2_ini+(b.X_co2_high-b.X_co2_low)*ss*randn(1,1);           
        end
        
        % compute predicted mass of co2 and h2p
        [eps_x0, ~,~, ~]= crystal_fraction_silicic(T_0,P_0, eps_g_cand,X_co2_cand,mH2O,mCO2,...
            rho_m0,rho_x0,rho_g0);
        
        eps_m0 = 1-eps_g_cand-eps_x0;
        mm_gas = X_co2_cand*param.mm_co2+(1-X_co2_cand)*param.mm_h2o;
        
        M_co2_predicted = mco2_dissolved_sat(X_co2_cand,P_0,T_0)*eps_m0*V_0*rho_m0+...
            eps_g_cand*rho_g0*V_0*X_co2_cand*param.mm_co2/mm_gas;
        
        M_h2o_predicted = meq_water(X_co2_cand,P_0,T_0)*eps_m0*V_0*rho_m0+...
            eps_g_cand*rho_g0*V_0*(1-X_co2_cand)*param.mm_h2o/mm_gas;
        
        
        % compute error from initial guess
        err_co2_cand = abs(M_co2_predicted-M_co2_0)/M_co2_0;
        err_h2o_cand = abs(M_h2o_predicted-M_h2o_0)/M_h2o_0;
        err_tot_cand= err_co2_cand+ err_h2o_cand;      
        
        % compare candidate to initial
        if err_tot_cand < err_tot_ini
            eps_g_exp_ini=eps_g_exp_cand;
            eps_g_ini=10^eps_g_exp_ini;
            X_co2_ini = X_co2_cand;
            err_tot_ini=err_tot_cand;
        end
        
        eps_g_exp_trace(count) = eps_g_exp_ini;
        eps_g_trace(count) = eps_g_ini;
        X_co2_trace(count) = X_co2_ini;
        err_tot_trace(count) = err_tot_ini;
           
    end
    
    if count==maxcount_MT && err_tot_ini > err_tol
        disp('timed out')
        phase = 2;
        X_co20=0;
        eps_g0=0;
        mco2_diss=Conc_co2;
        [eps_x0, ~,~, ~]= crystal_fraction_silicic(T_0,P_0, 0,0,mH2O,mCO2,...
            rho_m0,rho_x0,rho_g0); % Initial crystal fraction assume zero gas
        
    else
        phase = 3;
        X_co20=X_co2_ini;
        eps_g0=eps_g_ini;
        [eps_x0, ~,~, ~]= crystal_fraction_silicic(T_0,P_0, eps_g0,X_co20,mH2O,mCO2,...
            rho_m0,rho_x0,rho_g0); % Initial crystal fraction assume zero gas
        
        mco2_diss = mco2_dissolved_sat(X_co20,P_0,T_0);
        
    end
    
end

% Function to compute mass fraction of co2 dissolved in melt at saturation
    function sol=mco2_dissolved_sat(X,P,T)
        
        P_MPa=P/1e6;
             
        % function & coefficients from Liu et al 2005
        Pc=P_MPa*X;
        Pw=P_MPa*(1-X);
        
        c1=5668;
        c2=-55.99;
        c3=0.4133;
        c4=2.041e-3;
        
        sol       = Pc*(c1+c2*Pw)/T+Pc*(c3*Pw^0.5+c4*Pw^1.5);
          
        sol = sol/1e6;
        
    end

    function sol=meq_water(X,P,T)
        P=P/1e6;
        Pw=(1-X)*P;
        Pc=X*P;
        T_C=T;
        a1=354.94;
        a2=9.623;
        a3=-1.5223;
        a4=1.2439e-3;
        a5=-1.084e-4;
        a6=-1.362e-5;
        sol=(a1*Pw^0.5+a2*Pw+a3*Pw^1.5)/T_C+a4*Pw^1.5+Pc*(a5*Pw^0.5+a6*Pw);
        sol=sol/100;
    end


end

