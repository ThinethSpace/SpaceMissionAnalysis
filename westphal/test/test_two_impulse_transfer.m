
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

allowed_angle = deg2rad(15);
num_node_crossings = 15;


%[x_opt_all, t_w_all, eta_all] = OM.perform_three_impulse_transfer(a_leo, i_leo, RAAN_leo, a_geo, 0, lambda_t0, u_t0, mu, allowed_angle, num_node_crossings);
%score = sqrt((x_opt_all/norm(x_opt_all)).^2 + (t_w_all/norm(t_w_all)).^2);

[x_all, t_w_all, delta_v_all, delta_lambda_all, nu_opt_all] = OM.perform_two_impulse_transfer(a_leo, i_leo, RAAN_leo, a_geo, 0, lambda_t0, u_t0, mu, 40);


figure
subplot(5,1,1), plot(x_all), title('x')
subplot(5,1,2), plot(t_w_all), title('t_w')
subplot(5,1,3), plot(delta_v_all), title('\Delta v')
subplot(5,1,4), plot(delta_lambda_all), title('\Delta \lambda')
subplot(5,1,5), plot(nu_opt_all), title('\nu')

%nu = linspace(2.4, pi, 100);
%t_trans = zeros(size(nu));  % preallocate
%
%for k = 1:length(nu)
%    [~, ~, ~, ~, t_trans(k),~] = OM.one_tangent_burn_LEO_GEO(a_leo, a_geo, nu(k), mu);
%end
%
%plot(nu, t_trans)
%xlabel('\nu (rad)')
%ylabel('\phi_{trans}')
%grid on

