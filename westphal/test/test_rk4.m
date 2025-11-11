%%%%% Test of RK4 method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Instantiation function collector 

fc = FunctionCollector;

% derivative function

mu   = 398600.4418;          % [km^3/s^2]
T = 5.8262e+03;
h = 60;
Re   = 6378.1363;            % [km]
J_2 = 0.00108263;


f = @(t, rr, zz) zz;
% g = @(t, rr, zz) -(mu/(norm(rr)^3)) * rr;

g = @(t, rr, zz) ...
    [
    -(mu*rr(1)) / (norm(rr)^3) * (1 - (3/2)*J_2 * ((Re/norm(rr))^2 * (5*((rr(3)^2)/norm(rr)^2) - 1)));
    -(mu*rr(2)) / (norm(rr)^3) * (1 - (3/2)*J_2 * ((Re/norm(rr))^2 * (5*((rr(3)^2)/norm(rr)^2) - 1)));
    -(mu*rr(3)) / (norm(rr)^3) * (1 - (3/2)*J_2 * ((Re/norm(rr))^2 * (5*((rr(3)^2)/norm(rr)^2) - 3)));
    ];

rr0 = 1E3 * [6.9281; 0; 0];
vv0= [0; -1.0431; 7.5512];

[tt, R, V] = fc.perform_rk4_2nd_ODE(f, g, 0, rr0, vv0, 30 * T, h);

figure('Name','3D Orbit (ECI)');
plot3(R(:,1), R(:,2), R(:,3), 'LineWidth', 1.0); hold on; grid on; axis equal;
% Draw a translucent Earth
[fX,fY,fZ] = sphere(80);
surf(Re*fX, Re*fY, Re*fZ, 'FaceAlpha', 0.1, 'EdgeColor', 'none'); colormap gray;
xlabel('x_{ECI} [km]'); ylabel('y_{ECI} [km]'); zlabel('z_{ECI} [km]');
title('Perturbed Keplerian Orbit (ECI/ICRF)');
view(35,25);