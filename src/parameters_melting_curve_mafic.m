function [a dadx  dady dadz b dbdx dbdy dbdz]=parameters_melting_curve_mafic(mH2O,mCO2,P)

x=mH2O;
y=mCO2;
z=P/1e6;

intercept=-0.0106007180771044;
H2Ocoeff=0.00642879427079997;
CO2coeff=-0.000362698886591479;
Pcoeff=6.33762356329763e-06;
H2OxCO2coeff=-0.0000409319695736593;
H2OxPcoeff=-0.0000020971242285322;
CO2xPcoeff=3.66354014084072e-07;
H2Osquarecoeff=-0.00127225661031757;
CO2squarecoeff=0.000219802992993448;
Psquarecoeff=-1.4241625041626E-09;

a=intercept + H2Ocoeff*x + CO2coeff*y + Pcoeff*z + H2OxCO2coeff*x*y + ...
    H2OxPcoeff*x*z + CO2xPcoeff*y*z + H2Osquarecoeff*x^2 + CO2squarecoeff*y^2 + Psquarecoeff*z^2;

dadx= H2Ocoeff+H2OxCO2coeff*y+H2OxPcoeff*z+2*H2Osquarecoeff*x;

dady= CO2coeff+H2OxCO2coeff*x+CO2xPcoeff*z+2*CO2squarecoeff*y;

dadz= Pcoeff+H2OxPcoeff*x+CO2xPcoeff*y+2*Psquarecoeff*z;

dadx=dadx*100;
dady=dady*100;
dadz=1e-6*dadz;

clear intercept
clear H2Ocoeff
clear CO2coeff
clear Pcoeff
clear H2OxCO2coeff
clear H2OxPcoeff
clear CO2xPcoeff
clear H2Osquarecoeff
clear CO2squarecoeff
clear Psquarecoeff

%b value
intercept=12.1982401917454;
H2Ocoeff=-7.49690626527448;
CO2coeff=0.398381500262876;
Pcoeff=-0.00632911929609247;
H2OxCO2coeff=0.0571369994114008;
H2OxPcoeff=0.00216190962922558;
CO2xPcoeff=-0.000409810092770206;
H2Osquarecoeff=1.48907741502382;
CO2squarecoeff=-0.251451720536687;
Psquarecoeff=1.36630369630388e-06;

b=intercept + H2Ocoeff*x + CO2coeff*y + Pcoeff*z + H2OxCO2coeff*x*y + H2OxPcoeff*x*z + CO2xPcoeff*y*z + H2Osquarecoeff*x^2 + CO2squarecoeff*y^2 + Psquarecoeff*z^2;

dbdx= H2Ocoeff+H2OxCO2coeff*y+H2OxPcoeff*z+2*H2Osquarecoeff*x;

dbdy= CO2coeff+H2OxCO2coeff*x+CO2xPcoeff*z+2*CO2squarecoeff*y;

dbdz= Pcoeff+H2OxPcoeff*x+CO2xPcoeff*y+2*Psquarecoeff*z;

dbdx=dbdx*100;
dbdy=dbdy*100;
dbdz=1e-6*dbdz;

