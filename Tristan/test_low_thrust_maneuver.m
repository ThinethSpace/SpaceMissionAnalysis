%{
Input:
- x0          initial state
- tw          waiting time
- c[1..N]     throttles
- delta[1..N] steering angles

1) Coast phase:
   propagate x0 → x1 over tw with ac = 0

2) Thrust phase:
   for k = 1..N
       propagate over Δt_arc
       using controls (ck, δk)

Output:
- final state x(tf)
%}
OP = OrbitPropagation;

% Convert Kepplerian initial orbit into cartesian
a = 200 + 6378;     % km
e = 0;              % circular orbit
i = deg2rad(22.5);  % rad
OM= deg2rad(40);    % rad
nu= deg2rad(17);    % rad  
om= deg2rad(0);     % rad

% Information about Target in GEO
lambda_t0 = deg2rad(200); % rad

% Use given initial guess
params.N     = 12;
params.mu    = 398600.4418;             % km^3/s^2
params.amax  = 1e-7;                    % km/s^2
params.tw    = 3 * 3600;                % 3 hours (seconds)
params.tLT   = 10 * 86400;              % 10 days (seconds)
params.dt    = params.tLT / params.N ;  % s

[ params.r0, params.v0 ] = OP.convert_kep2car(a, e, i, OM, om, nu, params.mu );

params.rGEO = 42164;
params.nGEO = sqrt(params.mu/params.rGEO^3);

params.lambda_t0 = lambda_t0;

params.wa = 1;
params.we = 1;
params.wi = 1;
params.wl = 10;
params.wu = 1e-3;

alpha0 =  zeros(params.N,1);         % only tangential burns
delta0 = -0.05 * ones(params.N,1);   % small normal bias (radians)
c0     =  0.8  * ones(params.N,1);   % moderate–high throttle

% Initial Z0
z0 = [delta0; c0]; % Ignore alpha0

lb = [ -2/pi * ones(params.N,1);   % beta lower boundary
        0    * ones(params.N,1) ]; % thrust lower boundary

ub = [  2/pi * ones(params.N,1);   % beta upper boundary
        1    * ones(params.N,1) ]; % thrust upper boundary

%% let fmincom do its magic
OP = OrbitPropagation;

costfun = @(z) OP.objective_low_thrust(z, params);

opts = optimoptions('fmincon','Display','iter','Algorithm','sqp', ...
                    'MaxFunctionEvaluations',1e5);

[z_opt, J_opt] = fmincon(costfun, z0, [], [], [], [], lb, ub, [], opts);


figure;
subplot(2,1,1);
stairs(1:params.N, rad2deg(z_opt.delta_opt), 'b', 'LineWidth', 1.5);
ylabel('\delta_k [deg]'); xlabel('Thrust arc k'); grid on;

subplot(2,1,2);
stairs(1:params.N, c_opt, 'r', 'LineWidth', 1.5);
ylabel('Throttle c_k'); xlabel('Thrust arc k'); grid on;