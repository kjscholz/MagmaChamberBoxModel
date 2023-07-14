function [eps_x, deps_x_dP, deps_x_dT,deps_x_deps_g, deps_x_dmco2_t, deps_x_dmh2o_t]= crystal_fraction_silicic(T,P, eps_g,X_co2,mH2O,mCO2,rho_m,rho_x,rho_g)

% convert T to Celcius from Kelvin
T=T-273;

%compute coefficients
[a dadx  dady dadz b dbdx dbdy dbdz c dcdx dcdy dcdz]=parameters_melting_curve_silicic(100*mH2O,100*mCO2,P);

if T < 400

    disp(T)
end    
    
    
eps_x=a*erfc(b*(T-c));

% derivatives
deps_x_deps_g=-1;
deps_x_dT=-2*a*b*exp(-b^2*(T-c)^2)/sqrt(pi); 
deps_x_dP=dadz*erfc(b*(T-c))-2*a*(T-c)/sqrt(pi)*exp(-b^2*(T-c)^2)*dbdz+2*a*b/sqrt(pi)*exp(-b^2*(T-c)^2)*dcdz;
deps_x_dmco2_t=dady*erfc(b*(T-c))-2*a*(T-c)/sqrt(pi)*exp(-b^2*(T-c)^2)*dbdy+2*a*b/sqrt(pi)*exp(-b^2*(T-c)^2)*dcdy;
deps_x_dmh2o_t=dadx*erfc(b*(T-c))-2*a*(T-c)/sqrt(pi)*exp(-b^2*(T-c)^2)*dbdx+2*a*b/sqrt(pi)*exp(-b^2*(T-c)^2)*dcdx;
