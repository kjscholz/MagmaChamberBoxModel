function dydz = odeChamber(time,y,        ...
    sw,...
    T_in,...
    param)

global storeTime storeTemp phase
global total_Mass M_h2o M_co2


P              = y(1);
T              = y(2);
eps_g          = y(3);
X_co2          = y(7);

%effective gas molar mass
m_g = param.mm_co2*X_co2+param.mm_h2o*(1-X_co2);

if ~isempty(storeTime)
    if storeTime(end)==time
        storeTime(end) = time;
        storeTemp(end) =  T;
    else
        storeTime = [storeTime; time];
        storeTemp = [storeTemp; T];
        
    end
else
    storeTime = [storeTime; time];
    storeTemp = [storeTemp; T];
    
end

cross=find(abs(diff(sign(diff(storeTime))))>0, 1); % number of time interutpions
% disp(storeTime)
% disp(cross)

if ~isempty(cross)
    cross_time=  storeTime(end);
    storeTime = [storeTime(storeTime<cross_time); cross_time];
    storeTemp = [storeTemp(storeTime<cross_time); storeTemp(end)];
end



V              = y(4);
dV_dP          = V/param.beta_r;
dV_dT          = -V*param.alpha_r;

rho_m          = y(5);
drho_m_dP      = rho_m/param.beta_m;
drho_m_dT      = -rho_m*param.alpha_m;

rho_x          = y(6);
drho_x_dP      = rho_x/param.beta_x;
drho_x_dT      = -rho_x*param.alpha_x;

[rho_g,drho_g_dP,drho_g_dT] = eos_g(P,T);

M_h2o    = y(9);%m_h2o_rep_start + Mdot_in_pass*tot_h2o_frac_in*(time-begin_time);
M_co2    = y(10);%m_co2_rep_start + Mdot_in_pass*tot_co2_frac_in*(time-begin_time);
total_Mass = y(8);%Mass_tot_rep_start +  Mdot_in_pass*(time-begin_time);

m_h2o= M_h2o/(total_Mass);
m_co2=M_co2/(total_Mass);
%display(['Line 84 odeChamber : M_h2o = ' num2str(M_h2o) ' ; M_co2 = ' num2str(M_co2) ' ; m_h2o =  ' num2str(m_h2o) ' ; m_co2 = ' num2str(m_co2)])

if strcmp(param.composition,'silicic')
    [eps_x, deps_x_dP, deps_x_dT,deps_x_deps_g, deps_x_dmco2_t, deps_x_dmh2o_t]= crystal_fraction_silicic(T,P,eps_g, X_co2,m_h2o,m_co2,rho_m,rho_x,rho_g);%,m_h2o,m_co2,rho_m,rho_x,rho_g
elseif strcmp(param.composition,'mafic')
    [eps_x, deps_x_dP, deps_x_dT,deps_x_deps_g, deps_x_dmco2_t, deps_x_dmh2o_t]= crystal_fraction_mafic(T,P,eps_g, X_co2,m_h2o,m_co2,rho_m,rho_x,rho_g);%,m_h2o,m_co2,rho_m,rho_x,rho_g
end

eps_m=1-eps_x-eps_g;

rho            = eps_m*rho_m + eps_g*rho_g + eps_x*rho_x;

    if strcmp(param.composition,'silicic')
        [m_eq,dm_eq_dP,dm_eq_dT,dm_eq_dX_co2,C_co2_t,dC_co2_dP, dC_co2_dT, dC_co2_dX_co2] = exsolve_silicic(P,T, X_co2);
    elseif strcmp(param.composition,'mafic')
        [m_eq,dm_eq_dP,dm_eq_dT,dm_eq_dX_co2,C_co2_t,dC_co2_dP, dC_co2_dT, dC_co2_dX_co2] = exsolve_mafic(P,T, X_co2);
    end
    
if phase == 3
       
    C_co2=C_co2_t;   
        
elseif phase == 2   
 
    C_co2 = m_co2;
    
end


[rho_g,drho_g_dP,drho_g_dT] = eos_g(P,T);

rho            = eps_m*rho_m + eps_g*rho_g + eps_x*rho_x;
drho_dP        = eps_m*drho_m_dP + eps_g*drho_g_dP + eps_x*drho_x_dP+(rho_x-rho_m)*deps_x_dP;
drho_dT        = eps_m*drho_m_dT + eps_g*drho_g_dT + eps_x*drho_x_dT+(rho_x-rho_m)*deps_x_dT;
drho_deps_g    = -rho_m + rho_g;
drho_dX_co2    = 0;%-rho_m + rho_x;
drho_dMco2t =  (rho_x-rho_m)*deps_x_dmco2_t;% COMPUTE deps_x_dMco2t
drho_dMh2ot =  (rho_x-rho_m)*deps_x_dmh2o_t;% COMPUTE deps_x_dMh2ot


% specific heat of gas
[c_g, dc_g_dX_co2] = gas_heat_capacity(X_co2);

c              = (1/rho)*(rho_x*eps_x*param.c_x+rho_m*eps_m*param.c_m+rho_g*eps_g*c_g);

% computing the product of density and specific heat for the mixture and
% its derivatives
rc              = rho_x*eps_x*param.c_x+rho_m*eps_m*param.c_m+rho_g*eps_g*c_g;
drc_dP          = eps_x*param.c_x*drho_x_dP+eps_g*c_g*drho_g_dP+eps_m*param.c_m*drho_m_dP+(rho_x*param.c_x-rho_m*param.c_m)*deps_x_dP;
drc_dT          = eps_x*param.c_x*drho_x_dT+eps_g*c_g*drho_g_dT+eps_m*param.c_m*drho_m_dT+(rho_x*param.c_x-rho_m*param.c_m)*deps_x_dT;
drc_deps_g      = rho_g*c_g-rho_m*param.c_m;
drc_dXco2       = 0;

% boundary conditions
[Mdot_in, Mdot_out, Mdot_v_in, Mdot_v_out,Mdot_c_in, Mdot_c_out, Hdot_in, Hdot_out, P_loss,eta_r] ...
    = boundary_conditions(P,T,eps_g,eps_x,V,rho_m,rho_x,m_eq,rho_g,rho,c,sw,...
    X_co2, C_co2,T_in, param);


% coefficients in the system of unknowns Ax = B, here x= [dP/dt dT/dt deps_g/dt dX_co2/dt]
% note: P, T, and phi are y(1), y(2) and y(3) respectively
% values matrix A

% conservation of (total) mass
a11 = (1/rho)*drho_dP +(1/V)*dV_dP;
a12 = (1/rho)*drho_dT     + (1/V)*dV_dT ;%+ (1/rho)*drho_deps_x*deps_x_dT;
a13 = (1/rho)*drho_deps_g;
a14 = 0;

% conservation of water mass
a21 = eps_m*dm_eq_dP-m_eq*deps_x_dP+m_eq*eps_m*dV_dP/V+m_eq*eps_m* drho_m_dP/rho_m+(1-X_co2)/m_g*eps_g*param.mm_h2o*drho_g_dP/rho_m+...
    (1-X_co2)/m_g*eps_g*rho_g*param.mm_h2o*dV_dP/(rho_m*V);
a22 = eps_m*dm_eq_dT-m_eq*deps_x_dT+m_eq*eps_m*(1/V)*dV_dT+m_eq*eps_m*drho_m_dT/rho_m+(1-X_co2)/m_g*eps_g*param.mm_h2o*drho_g_dT/rho_m+...
    (1-X_co2)/m_g*eps_g*rho_g*param.mm_h2o*dV_dT/(rho_m*V);
a23 = -m_eq+(1-X_co2)/m_g*rho_g*param.mm_h2o/rho_m;
a24 = eps_m*dm_eq_dX_co2-1/m_g*eps_g*rho_g*param.mm_h2o/rho_m-(1-X_co2)*eps_g*rho_g*param.mm_h2o*(param.mm_co2-param.mm_h2o)/(m_g^2*rho_m);


% conservation of (total) enthalpy - all divided by rc*T*V

a31  = drc_dP/(rc)+dV_dP/V+param.L_m*(-eps_x*drho_x_dP/(rc*T)-rho_x*deps_x_dP/(rc*T)-rho_x*eps_x*dV_dP/(rc*V*T))+...
    param.L_e*(-dm_eq_dP*rho_m*eps_m/(rc*T)-m_eq*eps_m*drho_m_dP/(rc*T)+m_eq*rho_m*deps_x_dP/(rc*T)-m_eq*rho_m*eps_m*dV_dP/(rc*V*T));
a32  = drc_dT/(rc)+1/T+dV_dT/V+param.L_m*(-eps_x*drho_x_dT/(rc*T)-rho_x*deps_x_dT/(rc*T)-rho_x*eps_x*dV_dT/(rc*T*V))+...
    param.L_e*(-dm_eq_dT*rho_m*eps_m/(rc*T)-m_eq*eps_m*drho_m_dT/(rc*T)+m_eq*rho_m*deps_x_dT/(rc*T)-m_eq*rho_m*eps_m*dV_dT/(rc*T*V));
a33  = (rho_g*c_g-rho_m*param.c_m)/rc +param.L_e*m_eq*rho_m/(rc*T);
a34  = -param.L_e*rho_m*eps_m*dm_eq_dX_co2/(rc*T);


%conservation of CO2 mass

a41 = eps_m*dC_co2_dP-C_co2*deps_x_dP+C_co2*eps_m*dV_dP/V+C_co2*eps_m* drho_m_dP/rho_m+X_co2/m_g*eps_g*param.mm_co2*drho_g_dP/rho_m+...
    X_co2/m_g*eps_g*rho_g*param.mm_co2*dV_dP/(rho_m*V);
a42 = eps_m*dC_co2_dT-C_co2*deps_x_dT+C_co2*eps_m*(1/V)*dV_dT+C_co2*eps_m*drho_m_dT/rho_m+X_co2/m_g*eps_g*param.mm_co2*drho_g_dT/rho_m+...
    X_co2/m_g*eps_g*rho_g*param.mm_co2*dV_dT/(rho_m*V);
a43 = -C_co2+X_co2/m_g*rho_g*param.mm_co2/rho_m;
a44 = eps_m*dC_co2_dX_co2+1/m_g*eps_g*rho_g*param.mm_co2/rho_m-X_co2*eps_g*rho_g*param.mm_co2*(param.mm_co2-param.mm_h2o)/(m_g^2*rho_m);
% a42 = (param.mm_co2/m_g)*X_co2*(-param.alpha_r+(1/rho_g)*drho_g_dT)...
%     +((C_co2*rho_m*eps_m)/(eps_g*rho_g))*(-param.alpha_r-param.alpha_m-(1/eps_m)*deps_x_dT+(1/C_co2)*(dC_co2_dT));
% a43 = (param.mm_co2*X_co2)/(m_g*eps_g)-((C_co2*rho_m)/(eps_g*rho_g))*(1+deps_x_deps_g);
% a44 = (param.mm_co2*param.mm_h2o*X_co2)/(m_g^2)-(param.mm_co2^2/m_g^2)*X_co2+ param.mm_co2/m_g ...
%     -((C_co2*rho_m)/(eps_g*rho_g))*deps_x_dX_co2+(rho_m*eps_m/(rho_g*eps_g))*dC_co2_dX_co2;

% values vector B
% conservation of (total) mass

dM_h2o_t_dt=1/(rho*V)*((Mdot_v_in-Mdot_v_out)-m_h2o*(Mdot_in-Mdot_out));
dM_co2_t_dt=1/(rho*V)*((Mdot_c_in-Mdot_c_out)-m_co2*(Mdot_in-Mdot_out));


b1  =  (Mdot_in - Mdot_out)/(rho*V) - P_loss-(rho_x-rho_m)/rho*deps_x_dmh2o_t*dM_h2o_t_dt-(rho_x-rho_m)/rho*deps_x_dmco2_t*dM_co2_t_dt;
% conservation of water mass
b2  = (Mdot_v_in - Mdot_v_out)/(rho_m*V)+m_eq*(deps_x_dmh2o_t*dM_h2o_t_dt+deps_x_dmco2_t*dM_co2_t_dt)-m_eq*eps_m*P_loss-(1-X_co2)/m_g*rho_g/rho_m*eps_g*param.mm_h2o*P_loss;
%(Mdot_v_in - Mdot_v_out)/(rho_g*eps_g*V) - P_loss*((param.mm_h2o*(1-X_co2))/m_g+(m_eq*rho_m*eps_m)/(rho_g*eps_g));
% conservation of (total) enthalpy
b3  =  (Hdot_in - Hdot_out)/(rc*T*V) - 1/(rc*V*T)*((rho_x*param.c_x-rho_m*param.c_m)*T*V*(deps_x_dmh2o_t*dM_h2o_t_dt+deps_x_dmco2_t*dM_co2_t_dt))+...
    param.L_m*rho_x*V/(rc*V*T)*(deps_x_dmh2o_t*dM_h2o_t_dt+deps_x_dmco2_t*dM_co2_t_dt)+param.L_e*m_eq*rho_m*V/(rc*V*T)*(deps_x_dmh2o_t*dM_h2o_t_dt+deps_x_dmco2_t*dM_co2_t_dt)+...
    P_loss*(-1+param.L_m*rho_x*eps_x*V/(rc*V*T)+param.L_e*m_eq*rho_m*eps_m*V/(rc*V*T));%- P_loss*(1-(param.L_m*rho_x*eps_x)/(rho*c*T)-(param.L_e*m_eq*rho_m*eps_m)/(rho*c*T));
% conservation of CO2 mass
b4  =  (Mdot_c_in - Mdot_c_out)/(rho_m*V)+C_co2*(deps_x_dmh2o_t*dM_h2o_t_dt+deps_x_dmco2_t*dM_co2_t_dt)-C_co2*eps_m*P_loss-(X_co2)/m_g*rho_g/rho_m*eps_g*param.mm_co2*P_loss;

if phase == 3
    % set up matrices to solve using Cramer's rule
    A          = [ a11 a12 a13 a14; a21 a22 a23 a24; a31 a32 a33 a34; a41 a42 a43 a44];
    A_P        = [ b1  a12 a13 a14; b2  a22 a23 a24; b3  a32 a33 a34; b4  a42 a43 a44];
    A_T        = [ a11 b1  a13 a14; a21 b2  a23 a24; a31 b3  a33 a34; a41 b4  a43 a44];
    A_eps_g    = [ a11 a12 b1  a14; a21 a22 b2  a24; a31 a32 b3  a34; a41 a42 b4  a44];
    A_X_co2    = [ a11 a12 a13 b1 ; a21 a22 a23 b2 ; a31 a32 a33 b3 ; a41 a42 a43 b4 ];
    
    
    det_A          = det(A);
    
    dP_dt          = det(A_P)/det_A;
    dT_dt          = det(A_T)/det_A;
    deps_g_dt      = det(A_eps_g)/det_A;
    dX_co2_dt      = det(A_X_co2)/det_A;
    dV_dt          = dV_dP*dP_dt + dV_dT*dT_dt + V*P_loss;
    drho_m_dt      = drho_m_dP*dP_dt + drho_m_dT*dT_dt;
    drho_x_dt      = drho_x_dP*dP_dt + drho_x_dT*dT_dt;
else
    A          = [ a11 a12; a31 a32];
    A_P        = [ b1  a12;  b3  a32 ];
    A_T        = [ a11 b1 ; a31 b3 ];
    
    det_A          = det(A);
    dP_dt          = det(A_P)/det_A;
    dT_dt          = det(A_T)/det_A;
    deps_g_dt      = 0;
    
    dV_dt          = dV_dP*dP_dt + dV_dT*dT_dt + V*P_loss;
    drho_m_dt      = drho_m_dP*dP_dt + drho_m_dT*dT_dt;
    drho_x_dt      = drho_x_dP*dP_dt + drho_x_dT*dT_dt;
    dX_co2_dt      = 0;
end

dydz           = zeros(10,1);    % column vector
dydz(1)        = dP_dt;
dydz(2)        = dT_dt;
dydz(3)        = deps_g_dt;
dydz(4)        = dV_dt;
dydz(5)        = drho_m_dt;
dydz(6)        = drho_x_dt;
dydz(7)        = dX_co2_dt;
dydz(8)        = Mdot_in - Mdot_out;
dydz(9)        =  Mdot_v_in - Mdot_v_out;
dydz(10)       = Mdot_c_in - Mdot_c_out;

global eta_r_vec eta_t_vec
eta_r_vec = [eta_r_vec eta_r];
eta_t_vec = [eta_t_vec time];