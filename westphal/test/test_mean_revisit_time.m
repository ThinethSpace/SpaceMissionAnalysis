% Define constants
a = 7700; %km
e = 0.01; 
OM = 0;
om = 0;
nu0 = 0;
mu = 398600.440; %km³/s²
R_E = 6378;
J_2 = 0.00108263;

% Create instance of function collector
OP = OrbitPropagation;
OA = OrbitAnalysis;

% Calculate inclination for SSO orbit
% https://en.wikipedia.org/wiki/Sun-synchronous_orbit
i = acos(-(a / 12352)^(7/3));

% Calculate orbital period
T = 2 * pi * sqrt(a^3 / mu);

%  Propagate orbit and transform to keplarian elements
t0 = 0;
t1 = 30 * T;
t_step = 60;

% Propagate orbit
[tt, R, V, nunu, OmegaOmega] = OP.propagate_orbit_keplar_newton(a, e, i, OM, om, nu0, mu, t0, t1, t_step, R_E, J_2, false);

% Test
start_time = datetime([2010 1 17 10 20 36]);
tt_seconds = seconds(tt);

times = start_time + tt_seconds;

lla = eci2lla(R*1000, datevec(times));

[lat_grid, lon_grid] = meshgrid(-80:1:80, -180:1:180);
%[lla_grid, pts] = OA.create_grid(600);
%lat_grid = lla_grid(:,1); lon_grid = lla_grid(:,2);

mrt = OA.get_mean_revisit_time0(lat_grid, lon_grid, tt, lla, t_step, 1000);

mrt = mrt ./ 3600;

Utils.plot_colormap_earth([-80 80], [-180 180], mrt, 'Mean Revisit Time (Hours)')

%Utils.plot_orbit_3D(pts, 1, 'test', 'x', 'y', 'z', "dots")

Utils.plot_ground_track(lla(:,1), lla(:,2))
