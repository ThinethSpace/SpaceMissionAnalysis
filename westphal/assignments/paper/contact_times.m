tic
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
t1 = 60 * 60 * 24 * 365; % 1 year in seconds
t_step = 60;
t_start = datetime([2026 3 1 00 00 00]);

% Altiude and inclination ranges
altitudes = linspace(180, 400, 5); % km
inclinations = deg2rad(linspace(0, 100, 5)); 

%% Ground station grid and cloud probability

gs_rf = [
  37.00   -25.14  275.00
 -31.05   116.19  262.94
  -31.05  116.19  252.26
 -35.78   -69.40 1550.00
   5.25   -52.80  -14.67
  67.86    20.97  400.68
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


optical_probability_rf =[
  1.0;
  1.0;
  1.0;
  1.0;
  1.0;
  1.0;
  1.0;
  1.0;
  ];

optical_probability_optocom =[
  0.0;
  0.0;
  0.0;
  0.0;
  0.0;
  0.0;
  0.0;
  0.0;
    ];
%% Setup up simulation

% Create instance of function collector
OP = OrbitPropagation;
OA = OrbitAnalysis;
data = PassDuration;

lat_grid = gs_rf(:,1); lon_grid = gs_rf(:,2);
num_pts = length(lat_grid);

MIN_ELEV = 20; % deg

cloud_probability = optical_probability_rf;
min_elevation = (MIN_ELEV/180) * pi; % minimum elevation

results_pass_sum = zeros(length(altitudes), length(inclinations));
results_pass_count = zeros(length(altitudes), length(inclinations));



h = waitbar(0,'Propagation in progress...');
for i = 1:length(altitudes)
    a = R_E + altitudes(i);
    for j = 1:length(inclinations)
        waitbar(((i - 1) * length(inclinations) + j) / (length(altitudes) * length(inclinations)),h);


        inclination = inclinations(j);
        
        % Pass func
        % Init array for storing data
        data.last_access = false(num_pts,1);
        data.last_aos_time = nan(num_pts,1);
        data.pass_sum = zeros(num_pts,1);
        data.pass_count = zeros(num_pts,1);
        data.fom = zeros(num_pts,3);

        func = @(t, nu, Omega, rr, vv) OA.get_mean_pass_duration_increment(t_start, t, rr, lat_grid, lon_grid, min_elevation, R_E, cloud_probability,data);

        OP.propagate_orbit_keplar_newton(a, e, inclination, OM, om, nu0, mu, t0, t1, t_step, R_E, J_2, false, false, func);

        % Sovle for mean pass duration
        results_pass_sum(i, j) = sum(data.pass_sum / (60 * 60));
        results_pass_count(i, j) = sum(data.pass_count);
    end
end

close(h);               

pass_sum_table = array2table(results_pass_sum, 'VariableNames', string(inclinations), 'RowNames', string(altitudes));
pass_count_table = array2table(results_pass_count, 'VariableNames', string(inclinations ), 'RowNames', string(altitudes));
writetable(pass_sum_table, 'pass_sum_table.csv', 'WriteRowNames', true);
writetable(pass_count_table, 'pass_count_table.csv', 'WriteRowNames', true);
toc