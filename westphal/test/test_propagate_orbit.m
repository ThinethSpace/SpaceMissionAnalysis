% Initial keplarian orbits
a = 7000;
e = 0.01;
i = 10 / 180 * pi;
OM = 0;
om = 0;
nu0 = 0;
mu = 398600.4418;
t0 = 0;
t1 = 50;
R_E = 6378;
J_2 = 0.00108263;
% Propagate orbit
nu = propagate_orbit(a, e, i, OM, om, nu0, mu, t0, t1, R_E, J_2);

