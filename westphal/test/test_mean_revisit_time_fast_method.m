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
data = SharedData;

% Calculate inclination for SSO orbit
% https://en.wikipedia.org/wiki/Sun-synchronous_orbit
i = acos(-(a / 12352)^(7/3));

% Calculate orbital period
T = 2 * pi * sqrt(a^3 / mu);

%  Propagate orbit and transform to keplarian elements
t0 = 0;
t1 = 10 * T;
t_step = 60;
t_start = datetime([2010 1 17 10 20 36]);

% Create grid
[lla_grid, pts] = OA.create_grid(80000);
lat_grid = lla_grid(:,1); lon_grid = lla_grid(:,2);
num_pts = length(lat_grid);

% Pass func
% Init array for storing data
data.last_access = false(num_pts,1);
data.last_rise_time = nan(num_pts,1);
data.revisit_sum = zeros(num_pts,1);
data.revisit_count = zeros(num_pts,1);

func = @(t, nu, Omega, rr, vv) OA.get_mean_revisit_time_increment(t_start, t, rr, lat_grid, lon_grid,1000, data);

[tt, R, V, nunu, OmegaOmega] = OP.propagate_orbit_keplar_newton(a, e, i, OM, om, nu0, mu, t0, t1, t_step, R_E, J_2, false, true, func);

start_time = datetime([2010 1 17 10 20 36]);
tt_seconds = seconds(tt);
times = start_time + tt_seconds;
lla = eci2lla(R*1000, datevec(times));

% Sovle for mean revisit time
mrt = data.revisit_sum ./ data.revisit_count;
mrt(data.revisit_count == 0) = NaN;         

%Utils.plot_colormap_earth([-80 80], [-180 180], mrt, 'Mean Revisito Time (Hours)')
Utils.plot_ground_track(lla(:,1), lla(:,2))
Utils.plot_heatmap_earth(lat_grid, lon_grid, mrt, 'test')

