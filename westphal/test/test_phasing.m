% Earth constants
mu = 398600.44;       %km/s
Re = 6371;  %km

OM = OrbitManeuver;

% Interceptor LEO
h_leo = 200;    %km
a_leo = h_leo + Re;
i_leo = deg2rad(22.5);
RAAN_leo = deg2rad(40);
u_t0 = deg2rad(17);

% Target GEO
h_geo = 35786;  %km
a_geo = Re + h_geo;
lambda_t0 = deg2rad(200);


[delta_v_phase, delta_v_trans1, delta_v_trans2] = OM.perform_noncoplanar_phasing_circular_orbits(a_leo, i_leo, RAAN_leo, a_geo, 0, lambda_t0, u_t0, mu, 0, false)