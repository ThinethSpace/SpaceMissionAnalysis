% Create all orbits for an earth mars transfer

IT = InterplanetaryTransfers;
OP = OrbitPropagation;

mu_sun = 398600.440; %km³/s²;
mu_sun_converted = 39.4784176; %AU³/year²
mu_earth_converted = 1.1857 * 10^(-4);  %AU³/year²
perigee_height = 1.34 * 10E-6; %AU
Re = 4.26 * 10E-5; %AU

rr_1 = [0.605774717586, -0.80374565571, 0.00010684714]';
rr_2 = [-0.01300489410, 1.575580535344, 0.033161980867]';

delta_theta = deg2rad(143.451);
dt = 203 / 365.24219;

factors = [1, 1.001];

[vv_1, vv_2, a] = IT.solve_lamberts_problem_secant(rr_1, rr_2, delta_theta, dt, mu_sun_converted, factors, 20, 1E-10);

vv_earth = [4.909398453, 3.759466359, 0.0001209199]';
vv_mars = [-4.920097603, 0.411804813, 0.129389380]';

% Find position and velocity vector at hyperbolic perigee for orbit propagation
h = cross(vv_1, vv_2) / norm(cross(vv_1, vv_2));
vv_inf_earth = vv_1 - vv_earth;
v_perigee = sqrt(norm(vv_inf_earth)^2 + (2 * mu_earth_converted) / (Re + perigee_height));


function v_rot = rotateVector(v, k, theta)
    % ROTATEVECTOR rotates a vector v around axis k by angle theta (radians)
    % v: 3x1 vector to rotate
    % k: 3x1 rotation axis (must be unit vector)
    % theta: rotation angle in radians
    %
    % Returns v_rot: rotated vector

    % Ensure k is unit vector
    k = k / norm(k);

    % Rodrigues' rotation formula
    v_rot = v*cos(theta) + cross(k, v)*sin(theta) + k*dot(k, v)*(1 - cos(theta));
end

% Calculate values for hyperbolic earth orbit
r_p = Re + perigee_height;
e_earth = 1 + (r_p * norm(vv_inf_earth)^2) / mu_earth_converted;

delta = asin(1/e_earth) * 2;

% Rotate vv_inf_earth to perigee to find direction of velocity vector
vv_perige_norm = rotateVector(vv_inf_earth, h, -1*delta/2) / norm(vv_inf_earth);
vv_perige = vv_perige_norm * v_perigee;

rr_perige_norm = cross(vv_perige_norm, h) / norm(cross(vv_perige_norm, h));
rr_perige = rr_perige_norm * r_p;

[a, e, i, RAAN, omega, nu] = OP.convert_car2kep(rr_perige, vv_perige, mu_earth_converted);

t1 = 10 * 60 ;
t_step = 1 * 60;

Re_km = 6371;
Re = Re_km;

a_km = a * 1.495978707e8;
a = a_km

[tt, R, V] = OP.propagate_orbit_keplar_newton(a, e, i, RAAN, omega, nu, 398600, 0, t1, t_step, 0 , 0, true);


figure('Name','3D Orbit (ECI)');
plot3(R(:,1), R(:,2), R(:,3), 'LineWidth', 1.0); hold on; grid on; axis equal;
% Draw a translucent Earth
[fX,fY,fZ] = sphere(80);
surf(Re*fX, Re*fY, Re*fZ, 'FaceAlpha', 0.1, 'EdgeColor', 'none'); colormap gray;
xlabel('x_{ECI} [km]'); ylabel('y_{ECI} [km]'); zlabel('z_{ECI} [km]');
title('Unperturbed Keplerian Orbit (ECI/ICRF)');
view(35,25);