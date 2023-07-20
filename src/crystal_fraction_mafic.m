function [eps_x, deps_x_dP, deps_x_dT,deps_x_deps_g, deps_x_dmco2_t, deps_x_dmh2o_t]= crystal_fraction_mafic(T,P, eps_g,X_co2,mH2O,mCO2,rho_m,rho_x,rho_g)
% convert to C
T=T-273;
% compute coefficients
[a dadx  dady dadz b dbdx dbdy dbdz]=parameters_melting_curve_mafic(100*mH2O,100*mCO2,P);

eps_x=a*T+b;

% derivatives
deps_x_deps_g=-1;%this shouldn't be here...

deps_x_dT=a; 
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


dC_co2dXco2 = b3+b5*T+b7*P+2*b9*X_co2;
dC_co2dXco2 = 1e-6.*dC_co2dXco2;


if eps_x < 0 || eps_x>1
    eps_x=0;
    deps_x_dX_co2=0;
    deps_x_deps_g=0;
    deps_x_dT=0;
end
