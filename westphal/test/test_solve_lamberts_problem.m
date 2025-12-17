% Inputs
delta_theta = deg2rad(75);
dt = 180;
mu = 398600.440; %km³/s²

mu_converted = 0.00029591220828559; %AU³/days²


rr_1 = [1 0 0]';
rr_2 = [0.187, 0.698, 0]';

factors = [1, 2];

IT = InterplanetaryTransfers;

[vv_1, vv_2, a] = IT.solve_lamberts_problem_secant(rr_1, rr_2, delta_theta, dt, mu_converted, factors, 20, 1E-8);

fac = 1731.456;

vv_1 = vv_1 * fac
vv_2 = vv_2 * fac
