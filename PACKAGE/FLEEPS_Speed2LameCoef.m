function [Lambda, mu] = FLEEPS_Speed2LameCoef(C_l, C_t, rho)

% speed m/s
% density kg/m^3

mu = rho * C_t ^ 2;

Lambda = rho * ( C_l ^ 2 - 2 * C_t ^ 2);