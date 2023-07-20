function [Mdot_in, Mdot_out, Mdot_v_in, Mdot_v_out, Mdot_c_in, Mdot_c_out, Hdot_in, Hdot_out, P_loss,eta_r] ...
    = boundary_conditions(P,T,eps_g,eps_x,V,rho_m,rho_x,m_eq,rho_g,rho,c,sw,...
    X_co2,C_co2,T_in,param)
% This function passes in the boundary conditions for the ODE solver
global Mdot_in_pass Mdot_out_pass P_lit tot_h2o_frac_in tot_co2_frac_in
global Q_out_old M_h2o M_co2 total_Mass 

% set inflow conditions
rho_m_in       = rho_m; 
rho_x_in       = rho_x; 
P_in           = P;
eps_g_in       = 0.0;
X_co2_in       = 0.0;

[rho_g_in,~,~]  = eos_g(P_in,T_in);

if strcmp(param.composition,'silicic')
[eps_x_in, ~,~,~] = crystal_fraction_silicic(T_in,P_lit,eps_g_in, X_co2_in,tot_h2o_frac_in,tot_co2_frac_in,rho_m_in,rho_x_in,rho_g_in);%crystal_fraction(T,P, eps_g,X_co2,mH2O,mCO2,rho_m,rho_x,rho_g)
elseif strcmp(param.composition,'mafic')
[eps_x_in, ~,~,~] = crystal_fraction_mafic(T_in,P_lit,eps_g_in, X_co2_in,tot_h2o_frac_in,tot_co2_frac_in,rho_m_in,rho_x_in,rho_g_in);%crystal_fraction(T,P, eps_g,X_co2,mH2O,mCO2,rho_m,rho_x,rho_g)
end


rho_in         = (1-eps_g_in-eps_x_in)*rho_m_in + eps_g_in*rho_g_in + eps_x_in*rho_x_in;
[c_g_in, ~]    = gas_heat_capacity(X_co2_in);
c_in           = ((1-eps_g_in-eps_x_in)*rho_m_in*param.c_m + eps_g_in*rho_g_in*c_g_in + eps_x_in*rho_x_in*param.c_x)/rho_in;%c;

Mdot_in        = Mdot_in_pass;
Mdot_v_in      = tot_h2o_frac_in*Mdot_in;
Hdot_in        = c_in*T_in*Mdot_in;
Mdot_c_in      = tot_co2_frac_in*Mdot_in;


    
% set outflow conditions
if sw.eruption ==0
    Mdot_out       = 0;
    Mdot_v_out    = 0;
    Mdot_c_out   = 0;
elseif sw.eruption==1
    Mdot_out       = Mdot_out_pass;
    Mdot_v_out     = M_h2o/total_Mass*Mdot_out_pass;
    Mdot_c_out     = M_co2/total_Mass*Mdot_out_pass;
else
    disp('eruption not specified')
end

a                    = (V/(4*pi/3))^(1/3);   % chamber radius (m)
cc                   = 10*a; %outer shell radius (m)
dr                   = 0.1*a;
if ~isreal(a)
    disp('here')
    a                    = real((V/(4*pi/3))^(1/3));   % chamber radius (m)
    cc                   = 10*a; %outer shell radius (m)
    dr                   = 0.1*a;
end


if sw.heat_cond==1
    % heat loss
    Q_out          = heat_conduction_chamber(param.maxn,a,cc,dr,param.kappa,param.rho_r,param.c_r,param.Tb); 

elseif sw.heat_cond==0
    Q_out          = 0;
else
    disp('heat_cond not specified')
end

if isnan(Q_out)
    Q_out=Q_out_old;
    disp('Q_out is NaN AM I HERE?')
else
    Q_out_old=Q_out;
end

Hdot_out       = c*T*Mdot_out + Q_out;

% viscous relaxation
    [quadpts,weights] = GLQ_points_weights_hard(param.GLQ_n);
    b     = a +cc;
    quadpts_r = (b-a)/2.*quadpts + (a+b)/2;
    
      Trt = heat_conduction_chamber_profile(param.maxn,a,cc,quadpts_r,param.kappa,param.Tb);
A = param.A; % material-dependent constant for viscosity law (Pa s)
B = param.B; % molar gas constant (J/mol/K)
G = param.G; % activation energy for creep (J/mol)
eta_rt     = A.*exp(G./B./Trt);
I          = (b-a)/2.*sum(weights.*(eta_rt./(quadpts_r).^4));
eta_r      = 3*a^3*I;

if sw.visc_relax ==1
    P_loss = (P - P_lit)/eta_r;
elseif sw.visc_relax==0
    P_loss = 0;
else
    disp('visc_relax not specified')
end
