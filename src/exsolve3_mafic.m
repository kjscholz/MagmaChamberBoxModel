function [C_co2, X_co2] = exsolve3_mafic(P,T,m_eq) 

%convert to MPa and Celsius
P = P/1e6;
T = T-273.15;
m_eq=m_eq*1e2;

%NEWTON RAPHSON SOLVE FOR X_co2
errorTol = 1e-10;
h(1)=2.99622526644026;
h(2)=0.00322422830627781;
h(3)=-9.1389095360385;
h(4)=0.0336065247530767;
h(5)=0.00747236662935722;
h(6)=-0.0000150329805347769;
h(7)=-0.01233608521548;
h(8)=- 4.14842647942619e-6;
h(9)=- 0.655454303068124;
h(10)=- 7.35270395041104e-6;

%Water Paritioning Function
water=@(p,t,x,c) h(1)+h(2)*t+h(3)*x+h(4)*p+h(5)*t*x+h(6)*t*p+h(7)*x*p+...
    h(8)*t^2+h(9)*x^2+h(10)*p^2-c;

%(h(1)*(p*(1-x))^0.5+h(2)*(p*(1-x))+h(3)*(p*(1-x))^1.5)/t+...
%               h(4)*(p*(1-x))^1.5+p*x*(h(5)*(p*(1-x))^0.5+h(6)*(p*(1-x)))-c;
%Derivative of Water wrt Xco2
 
dwater_dx=@(p,t,x) h(3)+h(5)*t+h(7)*p+2*h(9)*x;
%P,T, and inital guesses/values
Xc_initial=0.01;
Xc_guess=Xc_initial;
Xc_prev=0;
count=0;
W=m_eq;
while abs(Xc_prev-Xc_guess)>errorTol %&& isreal(Xc_guess)
    count=count+1;
    if mod(count,1e6)==0
       % disp(count)
       % disp(Xc_guess)
    end
    Xc_prev=Xc_guess;
  %  Xc_guess= Xc_prev+(carbonDioxide(P,T,Xc_prev,Cc)/dcarbonDioxide_dx(P,T,Xc_prev));
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

%count;

X_co2=Xc_guess;

% MT adding this b/c some super CO2-rich cases make Xc_guess greater than 1
if X_co2>1
    X_co2=1
end


%partial pressures of CO2 and Water

%P_MPa=P/1e6;
Pc = P*X_co2;
Pw = P*(1-X_co2);
T = T+273.15; % convert back to Kelvin because that's what the Liu needs

% function & coefficients from Liu et al 2005


% C_co2 = 1.30173241482961 - 0.00231445700714088*T_C - ...
%     0.00325283455213666*X_co2 - 0.00209392520550635*P_MPa + ...
%     0.000114506368088678*T_C*X_co2 + 0.0000021271419946826*T_C*P_MPa +...
%     0.0000635046957136512*X_co2*P_MPa + 9.89390930903556e-7*T_C^2 - ...
%     0.0532032643851973*X_co2^2 + 7.2250847129467e-7*P_MPa^2;

% a=[1.9299e-3 6.7898e-4 1 1];
% 
% C_co2=(a(1)*Pc^0.5+a(2)*Pc)*erf(a(3)*(T_C-a(4)));

% a=[-1.3978e-5 1.0245e-3 1 1];
% C_co2=(a(1)*Pc^1.5+a(2)*Pc)*erf(a(3)*(T_C-a(4)));

c1=5668;
c2=-55.99;
c3=0.4133;
c4=2.041e-3;

C_co2       = Pc*(c1+c2*Pw)/T+Pc*(c3*Pw^0.5+c4*Pw^1.5);

C_co2       = C_co2*1e-6;