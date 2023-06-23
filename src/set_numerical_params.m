% Set numerical parameters

% error tolerances used in ode method
rel_tol        = 1e-7;
abs_tol        = 1e-7;

% maximum simulation time in seconds
max_simulation_time = 1e16; 

% Stuff for Meredith's IC Finder code
global ss err_tol maxcount_MT random_guess eps_g_exp_ini_pass X_co2_ini_pass
ss = 0.001; % step size used for MCMC in IC Finder
err_tol = 0.001; % error tolerance used for IC Finder
maxcount_MT = 2e5; % max count for IC Finder
random_guess=0; % start with a random guess? If not, pass initial guesses below
eps_g_exp_ini_pass = -3; % initial guess for the exponent of gas fraction
X_co2_ini_pass = 0.8; % initial guess for X_co2



