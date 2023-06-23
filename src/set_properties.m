
% Set recharging magma temperature 
% if you want recharging magma temperature to be the same as
% the initial magma chamber temp, write T_R = 'same'
% otherwise enter a temperature in Kelvin
T_R = 'same';

% Initial crystal fraction 
ini_eps_x  = 0.02;

% chamber eruption parameters
global P_crit
P_crit            = 20e6;   % critical overpressure (Pa)   
mdot_out          = 1e5;    % outflow rate during eruption (kg/s)
max_num_eruptions = 1000;   % set maximumnumber of eruptions
    
% magma properties
param.rho_m0       = 2400;      % initial melt density (kg/m^3)
param.beta_m       = 1e10;      % bulk modulus melt (Pa)
param.beta_x       = 1e10;      % bulk moduluis crystals (Pa)
param.alpha_m      = 1e-5;      % thermal expansion melt (K^-1)
param.alpha_x      = 1e-5;      % thermal expansion crystals (K^-1)
param.L_e          = 610e3;     % latent heat of exsolution (J kg^-1) value used by Caricchi and Blundy (2015)
param.mm_co2       = 44.01e-3;  % molar mass of CO2 (kg/mol)
param.mm_h2o       = 18.02e-3;  % molar mass of H2O (kg/mol)

if strcmp(param.composition,'silicic')
    param.rho_x0       = 2600;      % initial crystal density (kg/m^3)
    param.c_m          = 1200;      % specific heat of melt (J kg^-1 K-1)
    param.c_x          = 1200;      % specific heat of crystals (J kg^-1 K-1)
    param.L_m          = 290e3;     % latent heat of melting (J kg^-1) value used by Caricchi and Blundy (2015)
elseif strcmp(param.composition,'mafic')
    param.rho_x0       = 2900;      % initial crystal density (kg/m^3)
    param.c_m          = 1142;      % specific heat of melt (J kg^-1 K-1)
    param.c_x          = 1160;      % specific heat of crystals (J kg^-1 K-1)
    param.L_m          = 470e3;     % latent heat of melting (J kg^-1) value used by Caricchi and Blundy (2015)
end

% crustal properties
param.alpha_r        = 1e-5;      % thermal expansion country rock (K^-1)
param.beta_r         = 1e10;      % bulk modulus country rock (Pa)
param.kappa          = 1e-6;      % thermal conductivity (W m^-1 K^-1)
param.rho_r          = 2700;      % density crust (kg/m3)
param.c_r            = 1200;      % heat capacity crust (J kg^-1 K-1)

% Rheology of the crust:
param.A = 4.25e7; % material-dependent constant for viscosity law (Pa s)
param.B = 8.31; % molar gas constant (J/mol/K)
param.G = 141e3; % activation energy for creep (J/mol)
param.maxn = 10000;
param.GLQ_n =64;

% Geothermal gradient and surface temperature
T_surface    = 0+273; % surface temperature (K)
T_gradient   = 32/1e3; % thermal gradient (K/m)

% gravity
grav_acc      = 9.81; % gravitational acceleration (m/s2)  