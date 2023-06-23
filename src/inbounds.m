function [value] = inbounds(eps_g_exp_cand,X_co2_cand,b)

value = 0;

if eps_g_exp_cand > b.eps_g_exp_high
    value = 1;
end

if eps_g_exp_cand < b.eps_g_exp_low
    value = 1;
end

if X_co2_cand > b.X_co2_high
    value = 1;
end

if X_co2_cand < b.X_co2_low
    value = 1;
end

