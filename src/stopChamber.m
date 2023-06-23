function [value,isterminal,direction] = stopChamber(y,eruption,P_lit,P_crit,param)

P     = y(1);
T     = y(2);
eps_g = y(3);
V     = y(4);
X_co2 = y(7);
rho_m = y(5);
rho_x = y(6);
tot_m = y(8);
tot_w = y(9);
tot_c = y(10);

P_0    = P_lit;

param.mm_co2= 44.01e-3;
param.mm_h2o= 18.02e-3;
[rho_g,~,~]  = eos_g(P,T);

m_h20=tot_w/tot_m;
m_co2=tot_c/tot_m;

if strcmp(param.composition,'silicic')
[eps_x, ~,~,~]= crystal_fraction_silicic(T,P,eps_g,X_co2,m_h20,m_co2,rho_m,rho_x,rho_g);%crystal_fraction(T,P, eps_g,X_co2,mH2O,mCO2,rho_m,rho_x,rho_g)
elseif strcmp(param.composition,'mafic')
[eps_x, ~,~,~]= crystal_fraction_mafic(T,P,eps_g,X_co2,m_h20,m_co2,rho_m,rho_x,rho_g);%crystal_fraction(T,P, eps_g,X_co2,mH2O,mCO2,rho_m,rho_x,rho_g)
end


%% M. Townsend's new stuff to deal w/ transition from saturated/undersaturated

eps_m0=1-eps_x;

m_h2o_melt = tot_w/(V*rho_m*eps_m0);
m_co2_melt = tot_c/(V*rho_m*eps_m0);


if strcmp(param.composition,'silicic')
m_eq_max = exsolve_silicic(P, T, 0);
elseif strcmp(param.composition,'mafic')
m_eq_max = exsolve_mafic(P, T, 0);
end

value6 = m_h2o_melt-m_eq_max;

if strcmp(param.composition,'silicic')
[C_co2_sat,~]= exsolve3_silicic(P,T, m_h2o_melt);
elseif strcmp(param.composition,'mafic')
[C_co2_sat,~]= exsolve3_mafic(P,T, m_h2o_melt);
end
value8 = m_co2_melt-C_co2_sat;
    
%%

value1 = eps_x; % don't remelt to zero eps_x

value2 = eps_x./(1-eps_g)-0.8; % gas percolation threshold?

if eruption ==0
    value3 = (P-P_0)-P_crit; % pressure reached the critical
else
    value3 = -P_crit;
end
    
if eruption==1
    value4 = P_0-P; % pressure reached lithostatic again
else
value4 = -P_crit;
end

value5 = eps_x-0.5; % crystallizes to 50%

value7=P-(P_0-P_crit); % underpressure

value      = [value1; value2; value3; value4; value5; value6; value7; value8];

isterminal = [1; 1; 1; 1; 1; 1; 1;1];   % Stop the integration
direction  = [0; 0; 0; 0; 0; 0; 0;0];   % detect all zeros (default)