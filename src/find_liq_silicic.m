function Tl=find_liq_silicic(water,co2,P,ini_eps_x)

[a dadx  dady dadz b dbdx dbdy dbdz c dcdx dcdy dcdz]=parameters_melting_curve_silicic(100*water,100*co2,P);

%eps_x=a*erfc(b*(T-c));

f=@(x) a*erfc(b*(x-c))-ini_eps_x;
x0=1000; % in celsius
Tl=fzero(f,x0);
Tl=Tl+273.15;