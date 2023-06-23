function [meq,dmeqdP,dmeqdT,dmeqdXco2, C_co2, dC_co2dP, dC_co2dT, dC_co2dXco2] = exsolve_mafic(P,T,X_co2) 
%FIXED AS OF 2/09/2020

% compute water solubility and co2 solubility
P = P/1e6;
T_C=T-273.15;

b1=2.99622526644026;
b2=0.00322422830627781;
b3=-9.1389095360385;
b4=0.0336065247530767;
b5=0.00747236662935722;
b6=-0.0000150329805347769;
b7=-0.01233608521548;
b8=- 4.14842647942619e-6;
b9=- 0.655454303068124;
b10=- 7.35270395041104e-6;

meq       =b1+b2*T_C+b3*X_co2+b4*P+b5*T_C*X_co2+b6*T_C*P+b7*X_co2*P+...
    b8*T_C^2+b9*X_co2^2+b10*P^2;
dmeqdT    = b2+b5*X_co2+b6*P+2*b8*T_C;
dmeqdP    = b4+b6*T_C+b7*X_co2+2*b10*P;
dmeqdXco2 = b3+b5*T_C+b7*P+2*b9*X_co2;

meq       = 1e-2.*meq;
dmeqdP    = 1e-8.*dmeqdP;
dmeqdT    = 1e-2.*dmeqdT;
dmeqdXco2 = 1e-2.*dmeqdXco2;

%% Liu et al 2005 - rhyolite
c1=5668;
c2=-55.99;
c3=0.4133;
c4=2.041e-3;

Pc       = P*X_co2;
dPcdP    = X_co2;
dPcdXco2 = P;

Pw       = P*(1-X_co2);
dPwdP    = 1-X_co2;
dPwdXco2 = -P;

C_co2       = Pc*(c1+c2*Pw)/T+Pc*(c3*Pw^0.5+c4*Pw^1.5);
dC_co2dT    = -1*Pc*(c1+c2*Pw)*T^(-2);
dC_co2dP    = T^(-1)*(dPcdP*(c1+c2*Pw)+Pc*(c2*dPwdP))+dPcdP*(c3*Pw^0.5+c4*Pw^1.5)+...
              Pc*(0.5*c3*Pw^(-0.5)*dPwdP+1.5*c4*Pw^0.5*dPwdP);
dC_co2dXco2 = T^(-1)*(dPcdXco2*(c1+c2*Pw)+Pc*(c2*dPwdXco2))+dPcdXco2*(c3*Pw^0.5+c4*Pw^1.5)+...
               Pc*(0.5*c3*Pw^(-0.5)*dPwdXco2+1.5*c4*Pw^0.5*dPwdXco2);


C_co2       = 1e-6.*C_co2;
dC_co2dP    = 1e-12.*dC_co2dP;
dC_co2dT    = 1e-6.*dC_co2dT;
dC_co2dXco2 = 1e-6.*dC_co2dXco2;
