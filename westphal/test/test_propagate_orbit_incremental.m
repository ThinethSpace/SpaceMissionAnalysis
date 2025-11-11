% Define constants
a = 7000; %km
e = 0.01; 
OM = 0;
om = 0;
nu0 = 0;
mu = 398600.440; %km³/s²
R_E = 6378;
J_2 = 0.00108263;

% Create instance of function collector
fc = FunctionCollector;

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
[tt, R, V, nunu, OmegaOmega] = fc.propagate_orbit_keplar_newton(a, e, i, OM, om, nu0, mu, t0, t1, t_step, R_E, J_2);

figure('Name','3D Orbit (ECI)');
plot3(R(:,1), R(:,2), R(:,3), 'LineWidth', 1.0); hold on; grid on; axis equal;
% Draw a translucent Earth
[fX,fY,fZ] = sphere(80);
surf(R_E*fX, R_E*fY, R_E*fZ, 'FaceAlpha', 0.1, 'EdgeColor', 'none'); colormap gray;
xlabel('x_{ECI} [km]'); ylabel('y_{ECI} [km]'); zlabel('z_{ECI} [km]');
title('Perturbed Keplerian Orbit (ECI/ICRF)');
view(35,25);


