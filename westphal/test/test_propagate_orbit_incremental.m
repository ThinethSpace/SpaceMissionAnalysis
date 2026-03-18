% Define constants
a = 6500; %km
e = 0.01; 
OM = 0;
om = 0;
nu0 = 0;
mu = 398600.440; %km³/s²
R_E = 6378;
J_2 = 0.00108263;

% Create instance of function collector
OP = OrbitPropagation;

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

%figure('Name','3D Orbit (ECI)');
%plot3(R(:,1), R(:,2), R(:,3), 'LineWidth', 1.0); hold on; grid on; axis equal;
%[fX,fY,fZ] = sphere(80);
%surf(R_E*fX, R_E*fY, R_E*fZ, 'FaceAlpha', 0.1, 'EdgeColor', 'none'); colormap gray;
%xlabel('x_{ECI} [km]'); ylabel('y_{ECI} [km]'); zlabel('z_{ECI} [km]');
%title('Perturbed Keplerian Orbit (ECI/ICRF)');
%view(35,25);

% Test
start_time = datetime([2010 1 17 10 20 36]);
tt_seconds = seconds(tt);

times = start_time + tt_seconds;

lla = eci2lla(R*1000, datevec(times));

%Utils.plot_ground_track(lla(:,1), lla(:,2))
Utils.plot_orbit_3D(R, R_E, 'test', 'x', 'y', 'z', "lines")

