function Tl=find_liq_mafic(water,co2,P,ini_eps_x)

[a dadx  dady dadz b dbdx dbdy dbdz]=parameters_melting_curve_mafic(100*water,100*co2,P);

%ini_eps_x=a*T+b;
Tl=(ini_eps_x-b)/a;
Tl=Tl+273.15;