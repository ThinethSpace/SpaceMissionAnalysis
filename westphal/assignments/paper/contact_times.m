%% Constants and parameters
% Constants
R_E = 6378; % km
J_2 = 0.00108263;
mu = 398600.440; %km³/s²

% Non changing initial orbital parameters
e = 0.01; 
OM = 0;
om = 0;
nu0 = 0;

% Simulation parameters
t0 = 0;
t1 = 113 * T;
t_step = 60;
t_start = datetime([2010 1 17 10 20 36]);

% Define minimum elevation for contact
min_elevation = 1/180 * pi; % 10 degree minimum elevation

% Altiude and inclination ranges
altitude = linspace(200, 450, 50); % km
inclination = linspace(0, 100, 20); % degrees

%% Ground station grid and priority

gs_rf = [
    -25.14   37.00   275.00;
    5.15    50.00   386.68;
    116.19  -31.05   262.94;
    116.19  -31.05   252.26;
    -69.40  -35.78  1550.00;
    -52.80    5.25   -14.67;
    20.97   67.86   400.68;
    20.96   67.86   402.17;
    -4.37   40.45   794.09;
];

gs_optocom = [
 22.62   37.84   440
 25.13   35.33   400
 -16.51  28.30  2400
 -2.36   37.09   489
 10.13   52.85    72
 11.64   48.08   549
 13.07   53.33    66
 11.28   48.08   615];

%% Setup up simulation

% Create instance of function collector
OP = OrbitPropagation;
OA = OrbitAnalysis;
data = PassDuration;

lat_grid = gs_rf(:,1); lon_grid = gs_rf(:,2);
num_pts = length(lat_grid);

results = zeros(length(altitude), length(inclination));

for idx = 1:length(altitude)
    a = R_E + altitude(idx);
    for jdx = 1:length(inclination)
        i = inclination(jdx);
        
        % Pass func
        % Init array for storing data
        data.last_access = false(num_pts,1);
        data.last_rise_time = nan(num_pts,1);
        data.pass_sum = zeros(num_pts,1);
        data.pass_count = zeros(num_pts,1);

        func = @(t, nu, Omega, rr, vv) OA.get_mean_pass_duration_increment(t_start, t, rr, lat_grid, lon_grid, min_elevation, R_E, data);

        [tt, R, V, nunu, OmegaOmega] = OP.propagate_orbit_keplar_newton(a, e, i, OM, om, nu0, mu, t0, t1, t_step, R_E, J_2, false, false, func);

        % Sovle for mean pass duration
        mpd = data.pass_sum ./ data.pass_count;
        mpd(data.pass_count == 0) = NaN;         
    end
end