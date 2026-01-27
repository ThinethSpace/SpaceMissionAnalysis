
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

% Use given initial guess of chaser
params.N     = 12;
params.mu    = 398600.4418;             % km^3/s^2
params.amax  = 1e-7;                    % km/s^2
%params.amax = 0;
params.tw    = 3 * 3600;                % 3 hours (seconds)
params.tLT   = 10 * 86400;              % 10 days (seconds)
params.dt    = 60; % s 
params.nu0   = nu;                      % rad
params.n0    = sqrt(a ^ 3 / params.mu);
[ params.r0, params.v0 ] = OP.convert_kep2car(a, e, i, OM, om, nu, params.mu );

% Target initial information 
params.rGEO = 42164;
params.nGEO = sqrt(params.mu/params.rGEO^3);
params.lambda_t0 = lambda_t0;
params.e_GEO = 0;
params.i_GEO = 0;
params.OM_GEO = 0;
params.om_GEO = 0;

% Penalty weights
params.wa = 1;
params.we = 1;
params.wi = 1;
params.wl = 10;
params.wu = 0;

% Given initial guess for thrust profile
alpha0 =  zeros(params.N,1);         % only tangential burns
delta0 = -0.05 * ones(params.N,1);   % small normal bias (radians)
c0     =  0.8  * ones(params.N,1);   % moderate–high throttle

% Initial Z0
z0 = [delta0; c0]; % Ignore alpha0

% Lower and upper bounds for thrust profile
lb = [ -2/pi * ones(params.N,1);   % beta lower boundary
        0    * ones(params.N,1) ]; % thrust lower boundary

ub = [  2/pi * ones(params.N,1);   % beta upper boundary
        1    * ones(params.N,1) ]; % thrust upper boundary

%% let fmincom do its magic

%% ================== OPTIMIZATION ==================
costfun = @(z) OP.objective_low_thrust(z, params);

opts = optimoptions('fmincon', ...
    'Display','iter', ...
    'MaxIterations', 100, ...
    'MaxFunctionEvaluations', 20000, ...
    'StepTolerance', 1e-6);

[z_opt, J_opt] = fmincon(costfun, z0, [], [], [], [], lb, ub, [], opts);


disp('Optimization finished');
disp(['Final cost J = ', num2str(J_opt)]);



 %% ================== RE-PROPAGATE WITH OPTIMAL CONTROLS ==================
% This recomputes the full trajectory (coast + thrust)
[R_opt, V_opt, T_opt] = OP.propagate_low_thrust_history(z_opt, params);


%% ================== SANITY CHECKS ==================
% This is a sanity check, but something is not 100% correct
a = -params.mu / (2*(0.5*dot(V_opt(end,:),V_opt(end,:)) - params.mu/norm(R_opt(end,:))));

% Final time
tf = T_opt(end);

% Target final position
nu_t = params.lambda_t0 + params.nGEO * tf;

[r_t_final, ~] = OrbitPropagation.convert_kep2car( ...
    params.rGEO, params.e_GEO, ...
    params.i_GEO, params.OM_GEO, ...
    params.om_GEO, nu_t, params.mu);

% Miss distance
miss_distance = norm(R_opt(end,:) - r_t_final');
disp(['Final miss distance [km]: ', num2str(miss_distance)]);

% Longitude error
lambda_s = atan2(R_opt(end,2), R_opt(end,1));
lambda_t = params.lambda_t0 + params.nGEO * tf;
dLambda = wrapTo2Pi(lambda_s - lambda_t);


disp(['Longitude error [deg]: ', num2str(rad2deg(dLambda))]);

%% ================== PLOT TRAJECTORIES ==================
OP.plot_low_thrust_history(R_opt, T_opt, params);



%% =================== PLOT THRUST PROFILE ===================
figure;

delta_opt = z_opt(1:12);
c_opt = z_opt(13:24);
N = length(delta_opt);
arc_idx = 1:N;
% Subplot 1: Cos and Sin of steering angles
subplot(2,1,1);
plot(arc_idx, cos(delta_opt), '-o', 'LineWidth', 2, 'DisplayName', 'cos(\delta_k)');
hold on;
plot(arc_idx, sin(delta_opt), '-s', 'LineWidth', 2, 'DisplayName', 'sin(\delta_k)');
grid on;
xlabel('Thrust Arc k');
ylabel('Angle Component');
title('Steering Angle Components per Arc');
legend;


% Subplot 2: Throttle (stair plot)
subplot(2,1,2);
plot(arc_idx, c_opt, 'LineWidth', 2);
grid on;
xlabel('Thrust Arc k');
ylabel('Throttle c_k');
title('Throttle per Arc');
ylim([0 1]); % Because c_k ∈ [0,1]




