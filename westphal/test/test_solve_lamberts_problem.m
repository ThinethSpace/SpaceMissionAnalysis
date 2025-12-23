% Inputs Case 1
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

% Inputs Case 2

delta_theta = deg2rad(143.451);
dt = 203 / 365.24219;
mu_converted = 4*pi^2;

rr_1 = [0.605774717586, -0.80374565571, 0.00010684714]';
rr_2 = [-0.01300489410, 1.575580535344, 0.033161980867]';

factors = [1, 1.001];

[vv_1, vv_2, a] = IT.solve_lamberts_problem_secant(rr_1, rr_2, delta_theta, dt, mu_converted, factors, 20, 1E-10)
