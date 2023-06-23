function [C_co2, X_co2] = exsolve3_silicic(P,T,m_eq) 

%Takes pressure, temperature, and amount of water to solve for the
%concentration of CO2 and X_CO2 (basically, goes the other direction
%compared to exsolve.m) using a Newton Raphson scheme

%convert to MPa and Celsius
P = P/1e6;
m_eq=m_eq*1e2;

%NEWTON RAPHSON SOLVE FOR X_co2
errorTol = 1e-10;
h = [354.94  9.623 -1.5223 0.0012439 -1.084e-4 -1.362e-5];

%Water Paritioning Function
water=@(p,t,x,c) (h(1)*(p*(1-x))^0.5+h(2)*(p*(1-x))+h(3)*(p*(1-x))^1.5)/t+...
                h(4)*(p*(1-x))^1.5+p*x*(h(5)*(p*(1-x))^0.5+h(6)*(p*(1-x)))-c;
%Derivative of Water wrt Xco2
dwater_dx=@(p,t,x) -p*((1/T)*(0.5*h(1)*(p*(1-x))^(-0.5)+h(2)+1.5*h(3)*(p*(1-x))^0.5)+1.5*h(4)*(p*(1-x))^0.5+(p*x)*(0.5*h(5)*(p*(1-x))^(-0.5)+h(6)))...
            +p*(h(5)*(p*(1-x))^0.5+h(6)*(p*(1-x)));
        
%P,T, and inital guesses/values
Xc_initial=0.01;
Xc_guess=Xc_initial;
Xc_prev=0;
count=0;
W=m_eq;
while abs(Xc_prev-Xc_guess)>errorTol 
    count=count+1;
    if mod(count,1e6)==0
        disp(count)
        disp(Xc_guess)
    end
    Xc_prev=Xc_guess;
    Xc_guess= Xc_prev-(water(P,T,Xc_prev,W)/dwater_dx(P,T,Xc_prev));   
end

while isreal(Xc_guess)==0 && Xc_initial<=1
    Xc_initial=Xc_initial+0.01;
    Xc_prev=0;
    while abs(Xc_prev-Xc_guess)>errorTol %&& isreal(Xc_guess)
    count=count+1;
        if mod(count,1e6)==0
            disp(count)
            disp(Xc_guess)
        end
        Xc_prev=Xc_guess;
        %  Xc_guess= Xc_prev+(carbonDioxide(P,T,Xc_prev,Cc)/dcarbonDioxide_dx(P,T,Xc_prev));
         Xc_guess= Xc_prev-(water(P,T,Xc_prev,W)/dwater_dx(P,T,Xc_prev));   
    end
end

X_co2=Xc_guess;

%partial pressures of CO2 and Water
Pc       = P*X_co2;
Pw       = P*(1-X_co2);

% function & coefficients from Liu et al 2005
c1=5668;
c2=-55.99;
c3=0.4133;
c4=2.041e-3;

C_co2       = Pc*(c1+c2*Pw)/T+Pc*(c3*Pw^0.5+c4*Pw^1.5);
C_co2       = C_co2*1e-6;