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
i = 0;

% Calculate orbital period
T = 2 * pi * sqrt(a^3 / mu);

%  Propagate orbit and transform to keplarian elements
t0 = 0;
t1 = 4.3   * T;
t_step = 60;
t_start = datetime([2010 1 17 10 20 36]);

min_elevation = (20/180) * pi; % minimum elevation

% Create grid
gs_rf = [
  37.00   -25.14  275.00
  50.00    5.15   386.68
 -31.05   116.19  262.94
  %-31.05  %116.19  252.26
 -35.78   -69.40 1550.00
   5.25   -52.80  -14.67
  % 67.86  % 20.97  400.68
  67.86    20.96  402.17
  40.45    -4.37  794.09
];

gs_optocom = [
  37.84  22.62  440
  35.33  25.13  400
  28.30  -16.5 2400
  37.09  -2.36  489
  52.85  10.13   72
  48.08  11.64  549
  53.33  13.07   66
  48.08  11.28  615
];


test_optical_availability =[
  1.0;
  1.0;
  1.0;
  1.0;
  1.0;
  1.0;
  1.0;
  1.0;
  ];

gs = gs_rf;

lat_grid = gs(:,1); lon_grid = gs(:,2);
num_pts = length(lat_grid);

% Pass func
% Init array for storing data
data.last_access = false(num_pts,1);
data.last_aos_time = nan(num_pts,1);
data.pass_sum = zeros(num_pts,1);
data.pass_count = zeros(num_pts,1);
data.fom = zeros(num_pts,3);

func = @(t, nu, Omega, rr, vv) OA.get_mean_pass_duration_increment(t_start, t, rr, lat_grid, lon_grid, min_elevation, R_E, test_optical_availability, data);

[tt, R, V, nunu, OmegaOmega] = OP.propagate_orbit_keplar_newton(a, e, i, OM, om, nu0, mu, t0, t1, t_step, R_E, J_2, false, true, func);

start_time = datetime([2010 1 17 10 20 36]);
tt_seconds = seconds(tt);
times = start_time + tt_seconds;
lla = eci2lla(R*1000, datevec(times));

% Sovle for mean pass duration
mpd = data.pass_sum ./ data.pass_count;

mpd(data.pass_count == 0) = NaN;         

Utils.plot_ground_track(lla(:,1), lla(:,2))

hold on
plot(gs(:,2), gs(:,1), 'ro', 'MarkerFaceColor','g')  % ground stations
for n = 1:numel(gs(:,1))
    text(gs(n,2),gs(n,1),num2str(n), "FontSize", 14)
end

h_apo = a*(1+e) - R_E;
h_per = a*(1-e) - R_E;

altitude = R_E - norm(R(end,:));

Utils.plot_ground_station_range(deg2rad(gs), h_apo, h_per, min_elevation, R_E, true)

%Utils.plot_heatmap_earth(lat_grid, lon_grid, mpd, 'test')

