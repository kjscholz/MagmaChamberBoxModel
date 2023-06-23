function [meq,dmeqdP,dmeqdT,dmeqdXco2, C_co2, dC_co2dP, dC_co2dT, dC_co2dXco2] = exsolve_silicic(P,T,X_co2) 
% This script uses Liu et al. (2005) to calculate the solubility of water
% given some P,T, and X_CO2 conditions

% %partial pressures of CO2 and Water
P = P/1e6; % Have to convert to MPa

Pc       = P*X_co2;
dPcdP    = X_co2;
dPcdXco2 = P;

Pw       = P*(1-X_co2);
dPwdP    = 1-X_co2;
dPwdXco2 = -P;

b1=354.94;
b2=9.623;
b3=-1.5223;
b4=1.2439e-3;
b5=-1.084e-4;
b6=-1.362e-5;

meq       =(b1*Pw^0.5+b2*Pw+b3*Pw^1.5)/T+b4*Pw^1.5+Pc*(b5*Pw^0.5+b6*Pw);

dmeqdT    = -1*(b1*Pw^0.5+b2*Pw+b3*Pw^1.5)*T^(-2);

dmeqdP    = dPwdP*((1/T)*(0.5*b1*Pw^(-0.5)+b2+1.5*b3*Pw^0.5)+1.5*b4*Pw^0.5+Pc*(0.5*b5*Pw^(-0.5)+b6))...
            +dPcdP*(b5*Pw^0.5+b6*Pw);

dmeqdXco2 = dPwdXco2*((1/T)*(0.5*b1*Pw^(-0.5)+b2+1.5*b3*Pw^0.5)+1.5*b4*Pw^0.5+Pc*(0.5*b5*Pw^(-0.5)+b6))...
            +dPcdXco2*(b5*Pw^0.5+b6*Pw);

% coefficients for CO2 partitioning
c1=5668;
c2=-55.99;
c3=0.4133;
c4=2.041e-3;

C_co2       = Pc*(c1+c2*Pw)/T+Pc*(c3*Pw^0.5+c4*Pw^1.5);
dC_co2dT    = -1*Pc*(c1+c2*Pw)*T^(-2);
dC_co2dP    = T^(-1)*(dPcdP*(c1+c2*Pw)+Pc*(c2*dPwdP))+dPcdP*(c3*Pw^0.5+c4*Pw^1.5)+...
              Pc*(0.5*c3*Pw^(-0.5)*dPwdP+1.5*c4*Pw^0.5*dPwdP);
dC_co2dXco2 = T^(-1)*(dPcdXco2*(c1+c2*Pw)+Pc*(c2*dPwdXco2))+dPcdXco2*(c3*Pw^0.5+c4*Pw^1.5)+...
               Pc*(0.5*c3*Pw^(-0.5)*dPwdXco2+1.5*c4*Pw^0.5*dPwdXco2);

meq       = 1e-2.*meq;
dmeqdP    = 1e-8.*dmeqdP;
dmeqdT    = 1e-2.*dmeqdT;
dmeqdXco2 = 1e-2.*dmeqdXco2;

C_co2       = 1e-6.*C_co2;
dC_co2dP    = 1e-12.*dC_co2dP;
dC_co2dT    = 1e-6.*dC_co2dT;
dC_co2dXco2 = 1e-6.*dC_co2dXco2;
