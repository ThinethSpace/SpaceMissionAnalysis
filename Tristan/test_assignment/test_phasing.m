%% Task 1: Analytical phasing between LEO and GEO satellites
% Given parameters
h_km            = 200; % km
h_geo_km        = 35786; % km
i_deg           = 22.5; % deg
RAAN_deg        = 40; % deg
ut0_deg         = 17; % deg
mu_earth        = 398600.44; % km^3/s^2
R_earth         = 6371; % km
lambda_tar_0    = 200; % deg

%% a) Noncoplanar phasing maneuver. Calculate deltaV and trasfer time

% Step1: burn to phasing orbit
a_leo = R_earth + h_km; % km
a_geo = R_earth + h_geo_km; % km
omega_geo = sqrt(mu_earth / (a_geo)^3); % rad/s
omega_leo = sqrt(mu_earth / (a_leo)^3); % rad/s

% The transfer time can be calculated from half of the semimajor axis
% of the Hohman transfer orbit (elliptical orbit between LEO and GEO)
a_trans = (a_leo + a_geo) / 2; % km
t_trans = pi * sqrt(a_trans^3 / mu_earth);  % s
alpha_lead = omega_geo * t_trans;   % rad

% Calculate dTheta angle between ut0 and till RAAN. Also calculate Time
dTheta = pi - deg2rad(ut0_deg); % rad
dT_node = dTheta / omega_leo;  % s

% We want to calculate the new lead angle after the time to node
lambda_tar_1 = deg2rad(lambda_tar_0) + omega_geo * dT_node; % rad
lambda_int_1 = deg2rad(RAAN_deg + 180); % TODO: rad. 180 because   we are node point RAAN TODO
theta_new = lambda_int_1 - lambda_tar_1; % rad
alpha_new = theta_new + pi; % rad

% Now we determine the phasing orbit to determine it's semi major axis
P_phase = (alpha_new - alpha_lead) / omega_geo; % s

% Calculate semi major axis of phasing orbit
k_int = 1; % Only one revolution till node
a_phase = (mu_earth * (P_phase / (k_int * 2 * pi)) ^ 2 ) ^ (1/3); % km

% Calculate deltaV at Leo to enter phasing orbit
v_leo = sqrt(mu_earth / a_leo); % km
v_geo = sqrt(mu_earth / a_geo); % km
% Use vis-viva equation to calculate velocities
v_peri_phase = sqrt(2 * mu_earth / a_leo - mu_earth / a_phase); % km/s
v_trans1 = sqrt(2 * mu_earth / a_leo - mu_earth / a_trans); % km/s
v_trans2 = sqrt(2 * mu_earth / a_geo - mu_earth / a_trans); % km/s   

% calculate deltaVs
deltaV1 = abs(v_peri_phase - v_leo); % km/s
deltaV2 = abs(v_trans1 - v_peri_phase); % km/s
deltaVi = sqrt(v_geo ^ 2 + v_trans2 ^ 2 - 2 * v_geo * v_trans2  ... 
               * cos(deg2rad(i_deg))); % km/s
deltaV_total = deltaV1 + deltaV2 + deltaVi; % km/s

% total time is time to node + time till phase + transfer time
t_total = dT_node + P_phase + t_trans; % s

disp(['DeltaV to phasing orbit: ', num2str(deltaV1), ' km/s']);
disp(['DeltaV to transfer orbit: ', num2str(deltaV2), ' km/s']);
disp(['DeltaV combined circularazation and inclination change: ', num2str(deltaVi), ' km/s']); 
disp(['Total deltaV for non-coplanar phasing maneuver: ', num2str(deltaV_total), ' km/s']);
disp(['Transfer time for non-coplanar phasing maneuver: ', num2str(t_total), ' s']);
disp(' ');

%% b) Two burn inclination change at GEO

% Instead of doing combined maneuver of circularization and inclination change,
% we cando two separate burns at GEO. One circularaization burn and one inclination
% change burn.

deltaV_circ = abs(v_geo - v_trans2); % km/s
deltaV_incl = 2 * v_geo * sin(deg2rad(i_deg)    / 2); % km/s
deltaV_total_sep = deltaV1 + deltaV2 + deltaV_circ + deltaV_incl;
disp(['DeltaV for circularization at GEO: ', num2str(deltaV_circ), ' km/s']);
disp(['DeltaV for inclination change at GEO: ', num2str(deltaV_incl),   ' km/s']);
disp(['Total deltaV for two-burn inclination change at GEO: ', num2str(deltaV_total_sep), ' km/s']);