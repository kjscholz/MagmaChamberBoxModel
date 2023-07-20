%% silicic
T_s=1080-273;
P=200e6;
mH2O_s=linspace(.02,.06,20);
mCO2_s=1000e-6;
for i = 1: length(mH2O_s)

[a dadx  dady dadz b dbdx dbdy dbdz c dcdx dcdy dcdz]=parameters_melting_curve_silicic(mH2O_s(i),mCO2_s,P);
deps_x_dT(i)=-2*a*b*exp(-b^2*(T_s-c)^2)/sqrt(pi);
deps_x_dTdx(i)=-2*a*b*exp(-b^2*(T_s-c)^2)/sqrt(pi)*(-b^2*(T_s-c)^2);
end
figure
plot(mH2O_s, deps_x_dT_s)

hold on

%% mafic

T_m=1110-273;
P=200e6;
mH2O_m=linspace(.001,.03,20);
mCO2_m=500e-6;
for i = 1: length(mH2O_m)
[a dadx(i)  dady dadz b dbdx dbdy dbdz]=parameters_melting_curve_mafic(100*mH2O_m(i),100*mCO2_m,P);



deps_x_dmco2_t_m=dady*T_m+dbdy;
deps_x_dmh2o_t_m=dadx*T_m+dbdx;
deps_x_dT_m(i)=a;
end
plot(mH2O_m, dadx)