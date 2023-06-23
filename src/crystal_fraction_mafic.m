function [eps_x, deps_x_dP, deps_x_dT,deps_x_deps_g, deps_x_dmco2_t, deps_x_dmh2o_t]= crystal_fraction_mafic(T,P, eps_g,X_co2,mH2O,mCO2,rho_m,rho_x,rho_g)

%% NEW VERSION WITH SAGE's PARAMETERIZATION
T=T-273;

[a dadx  dady dadz b dbdx dbdy dbdz]=parameters_melting_curve_mafic(100*mH2O,100*mCO2,P);

eps_x=a*T+b;
% derivatives
deps_x_deps_g=-1;
deps_x_dT=a; % 

deps_x_dmco2_t=dady*T+dbdy;
deps_x_dmh2o_t=dadx*T+dbdx;
rho=eps_x*rho_x+eps_g*rho_g+(1-eps_x-eps_g)*rho_m;
dmtot_co2_dC_Co2=rho_m/rho*(1-eps_x-eps_g);
deps_x_dP=dadz*T+dbdz;

% Solubility for the derivatives (Liu et al, 2005)

P = P/1e6;


Pc       = P*X_co2;
dPcdP    = X_co2;
dPcdXco2 = P;

Pw       = P*(1-X_co2);
dPwdP    = 1-X_co2;
dPwdXco2 = -P;

b3=-0.00325283;
b5=0.000114506;
b7=0.000063505;
b9=-0.05320327;

% TK=T+273;
dC_co2dXco2 = b3+b5*T+b7*P+2*b9*X_co2;
dC_co2dXco2 = 1e-6.*dC_co2dXco2;


if eps_x < 0 || eps_x>1
    eps_x=0;
    deps_x_dX_co2=0;
    deps_x_deps_g=0;
    deps_x_dT=0;
end

% b              = 0.5;
% T_sol          = 700 + 273;
% T_liq          = 950 + 273;
% 
% 
% phi_x          = 1 - ((T-T_sol)/(T_liq-T_sol))^b;
% dphi_x_dT      = - b*(T-T_sol)^(b-1)/(T_liq-T_sol)^b;
% if T<T_sol
%     phi_x=1;
%     dphi_x_dT=0;
% elseif T>T_liq
%     phi_x=0;
%     dphi_x_dT=0;
% end
% eps_x          = (1-eps_g)*phi_x;
% deps_x_dT      = (1-eps_g)*dphi_x_dT;
% deps_x_deps_g  = -phi_x;
% deps_x_dX_co2 = 0;