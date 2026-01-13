% Create all orbits for an earth mars transfer

IT = InterplanetaryTransfers;
OP = OrbitPropagation;

% parse data files
[earth_data] = IT.parse_horizon_file('data/earth_data.txt', 'year'); % AU
[mars_data]  = IT.parse_horizon_file('data/mars_data.txt', 'year');   % AU

% Gravitational parameters
mu_sun   = 4*pi^2;        % AU^3 / year^2
mu_earth = 1.1857e-4;     % AU^3 / year^2
mu_mars  = 3.227e-5;      % AU^3 / year^2

% Radii of Planets
Re = 4.26e-5;             % AU
Rm = 2.27e-5;             % AU

% Periapsis altitude (200 km)
hp = 1.34e-6;             % AU


% Lambert Solver parameters
lsp.mu              = mu_sun;
lsp.factors         = [1, 1.001];
lsp.max_iterations  = 20;
lsp.tolerance       = 1E-8;

% Minimum solution for departure
[vv_inf_A, ephem_depA, ephem_arrA, dThetaA, dT_A] = IT.minimum_vinf_vals(...
    lsp, earth_data, mars_data, true);

% Plot Departure hyperbolic orbit
rrEarth_A = ephem_depA(1,2).Position; % AU
rrMars_A = ephem_arrA(1,2).Position;  % AU
vvEarth_A = ephem_depA(1,3).Velocity; % AU / Year
vvMars_A = ephem_arrA(1,3).Velocity;  % AU / Year

% Solve Lambert's problem for minimum departure C3
[vv1_A, vv2_A, ~] = IT.solve_lamberts_problem_secant(rrEarth_A, rrMars_A, dThetaA, dT_A, mu_sun, ...
                                lsp.factors, lsp.max_iterations, lsp.tolerance);


% Find position and velocity vector at hyperbolic perigee for orbit
% propagation
rEarth_p        = Re + hp; 
crossA          = cross(vvEarth_A, vvMars_A);    
hA              = crossA / norm(crossA);        % normalize
vv_inf_earth    = vv1_A - vvEarth_A;      
v_EarthPerigee  = sqrt(norm(vv_inf_earth)^2 + (2 * mu_earth) / rEarth_p); % Magnitude NOT VECTOR
e_earth         = 1 + (rEarth_p * norm(vv_inf_earth)^2) / mu_earth;
deltaA          = asin(1 / e_earth) * 2;

% Rotate vv_inf_earth to get perigee tangential velocity. Position vector
% results from the cross product between impulse and velocity
vv_perigee_earth_norm   = IT.rotateVector(vv_inf_earth, hA, -deltaA / 2) / norm(vv_inf_earth);
vv_perigee_earth        = vv_perigee_earth_norm * v_EarthPerigee;
rr_perigee_earth_norm   = cross(vv_perigee_earth_norm, hA) / norm(cross(vv_perigee_earth_norm, hA));
rr_perigee_earth        = rr_perigee_earth_norm * rEarth_p;

% Extract Kepler Elements and convert from AU / year -> km / s
[aE, eE, iE, RAANE, omegaE, nuE] = OP.convert_car2kep(rr_perigee_earth, vv_perigee_earth, mu_earth);
% ---- Unit conversion: AU/year -> km/s ----
AU2km   = 1.495978707e8;                 % km
yr2s    = 365.24219*24*3600;             % s
% Lengths Earth 
aE       = aE * AU2km;                   % km
Re       = Re * AU2km;                   % km
rr_perigee_earth = rr_perigee_earth * AU2km;           % km
% Velocities Earth 
vv_perigee_earh = vv_perigee_earth * AU2km / yr2s;    % km/s
% Time
t1 = 90 * 60; % s
t_step = 60;  % s
mu_earth_converted = 398600;    % km^3/s^2

% Propagate and Plot Earth hyperbolic departure
[tt, R, V, nunu, OMOM] = OP.propagate_orbit_keplar_newton(aE, eE, iE, RAANE, omegaE, nuE, ...
                            mu_earth_converted, 0, t1, t_step, 0 , 0, true);


figure('Name','Minimum C3 Departure 3D Orbit (ECI)');
plot3(R(:,1), R(:,2), R(:,3), 'LineWidth', 1.0); hold on; grid on; axis equal;
% Draw a translucent Earth
[fX,fY,fZ] = sphere(80);
surf(Re*fX, Re*fY, Re*fZ, 'FaceAlpha', 0.1, 'EdgeColor', 'none'); colormap gray;
xlabel('x_{ECI} [km]'); ylabel('y_{ECI} [km]'); zlabel('z_{ECI} [km]');
title('Departure Hyperbolic Orbit from Earth(ECI/ICRF)');
view(35,25);
