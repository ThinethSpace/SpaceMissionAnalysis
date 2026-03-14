% Define constants
e = 0.01; 
OM = 0;
om = 0;
nu0 = 0;
mu = 398600.440; %km³/s²
R_E = 6378;
J_2 = 0.00108263;

altitude = 2000; % km
a = R_E + altitude; % km

% Create instance of function collector
OP = OrbitPropagation;
OA = OrbitAnalysis;
data = PassDuration;

% Calculate inclination for SSO orbit
% https://en.wikipedia.org/wiki/Sun-synchronous_orbit
%i = acos(-(a / 12352)^(7/3));
i = pi/4;

% Calculate orbital period
T = 2 * pi * sqrt(a^3 / mu);

%  Propagate orbit and transform to keplarian elements
t0 = 0;
t1 = 10 * T;
t_step = 60;
t_start = datetime([2010 1 17 10 20 36]);

min_elevation = 1/180 * pi; % 10 degree minimum elevation

% Create grid
gs_rf = [
 -25.14   37.00   275.00
  5.15    50.00   386.68
 116.19  -31.05   262.94
 116.19  -31.05   252.26
 -69.40  -35.78  1550.00
 -52.80    5.25   -14.67
  20.97   67.86   400.68
  20.96   67.86   402.17
  -4.37   40.45   794.09
];

gs_optocom = [
 22.62   37.84   440
 25.13   35.33   400
 -16.51  28.30  2400
 -2.36   37.09   489
 10.13   52.85    72
 11.64   48.08   549
 13.07   53.33    66
 11.28   48.08   615
];

lat_grid = gs_rf(:,1); lon_grid = gs_rf(:,2);
num_pts = length(lat_grid);
priority = [1, 2, 3, 4]; % Priority for each ground station (1 = highest)

% Pass func
% Init array for storing data
data.last_access = false(num_pts,1);
data.last_rise_time = nan(num_pts,1);
data.pass_sum = zeros(num_pts,1);
data.pass_count = zeros(num_pts,1);

func = @(t, nu, Omega, rr, vv) OA.get_mean_pass_duration_increment(t_start, t, rr, lat_grid, lon_grid, min_elevation, R_E, data);

[tt, R, V, nunu, OmegaOmega] = OP.propagate_orbit_keplar_newton(a, e, i, OM, om, nu0, mu, t0, t1, t_step, R_E, J_2, false, true, func);

start_time = datetime([2010 1 17 10 20 36]);
tt_seconds = seconds(tt);
times = start_time + tt_seconds;
lla = eci2lla(R*1000, datevec(times));

% Sovle for mean pass duration
mpd = data.pass_sum ./ data.pass_count;
mpd(data.pass_count == 0) = NaN;         

Utils.plot_ground_track(lla(:,1), lla(:,2))
%Utils.plot_heatmap_earth(lat_grid, lon_grid, mpd, 'test')

