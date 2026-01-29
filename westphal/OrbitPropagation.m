classdef OrbitPropagation
    %FUNCTIONS Summary of this class goes here
    %   Detailed explanation goes here
    methods (Static)
        function [ rr, vv ] = convert_kep2car(a, e, i, OM, om, nu, mu )

            %%%%%%%Author: Elizabeth Klioner, TUB 2022, ALL RIGHTS RESERVED%%%%%%%%%%%%
            
            %%%% Input 
            
            % a [1x1] semi-major axis [km]
            % e [1x1] eccentricity [-]
            % i [1x1] inclination [rad]
            % OM [1x1] RAAN [rad]
            % om [1x1] argument of periapsis [rad]
            % nu [1x1] true anomaly [rad]
            % mu [1x1] gravitational parameter [km^3/s^2]
            
            %%%% Output
            
            % rr [3x1] position vector [km]
            % vv [3x1] velocity vector [km/s]
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            p = a * (1 - e^2); % semilatus rectum
            r = p / (1 + e*cos(nu)); % Radius 
            
            % Perifocal plane
            rr_pf = [ cos(nu); sin(nu); 0 ]*r; % radius in perifocal frame
            vv_pf = sqrt(mu/p)*[ -sin(nu); (e + cos(nu)); 0 ]; % velocity in perofical frame
            
            % Rotation matrices
            R_OM = [ cos(OM) sin(OM) 0; -sin(OM) cos(OM) 0; 0 0 1 ]; % Rotation around Om
            R_i = [ 1 0 0;  0 cos(i) sin(i); 0 -sin(i) cos(i) ]; % Rotation around i
            R_om = [ cos(om) sin(om) 0; -sin(om) cos(om) 0; 0 0 1 ]; % Rotation around om
            
            R = R_om*R_i*R_OM; % Total rotation matrix
            
            rr = R'*rr_pf; % position vector 
            vv = R'*vv_pf; % velocity vector
        
        end
        
        function [a, e, i, RAAN, omega, nu] = convert_car2kep(rr, vv, mu)
            rnorm = norm(rr);
            vnorm = norm(vv);

            h = cross(rr, vv);
            hnorm = norm(h);

            kvec = [0; 0; 1];
            i = acos(h(3) / hnorm);

            n = cross(kvec, h);
            nnorm = norm(n);

            evec = (1/mu) * ( (vnorm^2 - mu/rnorm) * rr - dot(rr, vv) * vv );
            e = norm(evec);

            RAAN = acos(n(1) / nnorm);
            if n(2) < 0
                RAAN = 2*pi - RAAN;
            end

            omega = acos(dot(n, evec) / (nnorm * e));
            if evec(3) < 0
                omega = 2*pi - omega;
            end

            nu = acos(dot(evec, rr) / (e * rnorm));
            if dot(rr, vv) < 0
                nu = 2*pi - nu;
            end

            a = 1 / (2/rnorm - vnorm^2/mu);
        end
        
        function [lat_deg, lon_deg] = convert_eci2lla(R_eci, tt, we, theta_g0)
            %%%%%%%Author: Kolja Westphal, TUB 2025, ALL RIGHTS RESERVED%%%%%%%%%%%%
            
            %%%% Input
            
            %%%% Output
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            N  = size(R_eci,1);
            lat_deg = zeros(N,1);
            lon_deg = zeros(N,1);
            
            t0 = tt(1);

            for k = 1:N
                theta = theta_g0 + we*(tt(k) - t0);   % Earth rotation since t0
                % ECI -> ECEF: rotate about z by -theta
                R3m = [ cos(theta)  sin(theta) 0;
                    -sin(theta)  cos(theta) 0;
                    0           0     1];
                r_ecef = R3m * R_eci(k,:).';
                
                x = r_ecef(1); y = r_ecef(2); z = r_ecef(3);
                r_norm = hypot(hypot(x,y), z);
                
                lat = asin(z / r_norm);
                lon = atan2(y, x);
                
                lat_deg(k) = rad2deg(lat);
                lon_deg(k) = rad2deg(lon);
            end
        end

        function [x, x_all] = perform_newtons_method(f, df, x0, opts)
            
            %%%%%%%Author: Kolja Westphal, ALL RIGHTS RESERVED%%%%%%%%%%%%
            
            %%%% Input 
            
            % f func function
            % df func derivate of function
            % x0 [1x1] start value for iteration [-]
            % opts optional struct optional arguments: 
            %   opts.k_max : maximum iterations
            %   opts.tol : stop criteria is abs(f(x)) < tol
            %   opts.return_all : returns estimates at all iterations if set to "true"
            
            %%%% Output
            
            % x [1x1] solution of x for f(x) = 0 [-]
            
            %%%% Example how to use:
            %  f = @(E) E - 0.5*sin(E) - 0.6
            % df = 1 - 0.5*cos(E)
            % sol = perform_newtons_method(f, df, 0.1)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % ----------------------------------
            % Sets (or defaults) solver options.
            % ----------------------------------

            % sets maximum number of iterations (defaults to 200)
            if (nargin < 4) || isempty(opts) || ~isfield(opts,'k_max')
                k_max = 200;
            else
                k_max = opts.k_max;
            end

            % sets tolerance (defaults to 10⁻¹⁰)
            if (nargin < 4) || isempty(opts) || ~isfield(opts,'tol')
                tol = 1e-12;
            else
                tol = opts.tol;
            end

            % determines if all intermediate estimates should be returned
            if (nargin < 4) || isempty(opts) || ~isfield(opts,'return_all')
                return_all = false;
            else
                return_all = opts.return_all;
            end



            % ----------------
            % Newton's method.
            % ----------------

            % returns initial guess if it is a root of f(x)
            if f(x0) == 0
                x = x0;
                return
            end

            % preallocates array
            if return_all
                x_all = zeros(1,k_max+1);
            end

            % iteration
            % Initialize the current guess for the iteration
            x_current = x0;

            for k = 1:k_max

                % stores results in arrays
                if return_all
                    x_all(k) = x_current;
                end

                f_current = f(x_current);
                if abs(f_current) < tol
                    x = x_current;
                    return
                end

                df_current = df(x_current);
                if df_current == 0
                    error('Derivative is zero. No solution found.');
                end

                % Update the current guess using Newton's method formula
                x_next = x_current - ( f_current / df_current );

                % Update current guess for the next iteration
                x_current = x_next;
            end

            % If maximum iterations reached without convergence
            warning('Maximum iterations reached. No solution found.');
            x = x_current;

        end
    
        function [xx, Y, Z] = perform_rk4_2nd_ODE(f, g, xx0, yy0, zz0, x1, h)
            
            %%%%%%%Author: Kolja Westphal, ALL RIGHTS RESERVED%%%%%%%%%%%%
            
            % Vector form of RK4 method
            
            %%%% Input 
            
            % f func derivative function for first first derivative
            % g func derivative function for second derivative
            % x0 [1x1] initial value for independent variable
            % y0 [3x1] initial vector for dependent vector
            % z0 [3x1] initial vector for first derivative of dep. vector
            % x1 [1x1] final value for dependent variable
            % h  [1x1] step size

            %%%% Output
            
            % Y, Z [nx3] matrix for extrapolated vectors
           
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Calculate steps
            Nx = floor((x1 - xx0)/h) + 1;
            xx = linspace(xx0, x1, Nx).';

            % Initial arrays
            Y = zeros(Nx, 3);
            Z = zeros(Nx, 3);

            % Calculate parameters and solve 
            for n = 1:Nx

                x_current = xx(n);

                k_0 = h * f(x_current, yy0, zz0);
                l_0 = h * g(x_current, yy0, zz0);
                k_1 = h * f(x_current + (1/2)*h, yy0 + (1/2)*k_0, zz0 + (1/2)*l_0);
                l_1 = h * g(x_current + (1/2)*h, yy0 + (1/2)*k_0, zz0 + (1/2)*l_0);
                k_2 = h * f(x_current + (1/2)*h, yy0 + (1/2)*k_1, zz0 + (1/2)*l_1);
                l_2 = h * g(x_current + (1/2)*h, yy0 + (1/2)*k_1, zz0 + (1/2)*l_1);
                k_3 = h * f(x_current + h, yy0 + k_2, zz0 + l_2);
                l_3 = h * g(x_current + h, yy0 + k_2, zz0 + l_2);
    
                yy1 = yy0 + (1/6)*(k_0 + 2*k_1 + 2*k_2 + k_3);
                zz1 = zz0 + (1/6)*(l_0 + 2*l_1 + 2*l_2 + l_3);
                
                Y(n,:) = yy1';
                Z(n,:) = zz1';

                % Set y and z for next iteration
                yy0 = yy1;
                zz0 = zz1;

            end

        end
    end
    methods
        function [nu1, OM1, rr, vv] = propagate_orbit_increment_keplar_newton(obj, a, e, i, Omega0, omega0, nu0, mu , t0, t1, R_E, J_2)

            %%%%%%%Author: Kolja Westphal, ALL RIGHTS RESERVED%%%%%%%%%%%%
            
            %%%% Input 
            
            % a [1x1] semi-major axis [km]
            % e [1x1] eccentricity [-]
            % i [1x1] inclination [rad]
            % Omega0 [1x1] RAAN [rad]
            % omega0 [1x1] argument of periapsis [rad]
            % nu0 [1x1] true anomaly [rad]
            % mu [1x1] gravitational parameter [km^3/s^2]
            % t0 [1x1] initial time [s]
            % t1 [1x1] time of progation [s]
            % R_E [1x1] body's equatorial radius [km]
            % J_2 [1x1] body's second dynamic form factor [-]
            
            %%%% Output
            
            % nu1 [1x1] propagated true anomaly
            % OM1 [1x1] propagated RAAN with J_2
            % rr [3x1] propagated position vector
            % vv [3x1] propagated velocity vector
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
            % Calculate initial values
        
            E0 = 2*atan2( sqrt(1-e)*sin(nu0/2), sqrt(1+e)*cos(nu0/2) );
        
            M0 = E0 - e * sin(E0);
            n = sqrt(mu/(a^3));
            M = M0 + n * (t1 - t0);
        
            % Define the function for Newton's method
            f = @(E) E - e*sin(E) - M;
            df = @(E) 1 - e*cos(E);
        
            % Initial guess for E
            if e > 0.8
                E_init = pi;
            else
                E_init = M;
            end
        
            % Run Newton's method and calcualte true anomaly
        
            E = wrapTo2Pi(obj.perform_newtons_method(f, df, E_init));
        
            nu1 = 2*atan2( sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2) );
            
            % Calculate node drift and transfrom to Cartesian Coordinates
        
            OM_dot = (-3/2) * sqrt(mu/(a^3)) * J_2 * (R_E / a)^2 * (cos(i)/(1-e^2)^2);
            OM1 = Omega0 + OM_dot * (t1 - t0);

            om_dot = (3/4)*J_2*n * (R_E/a)^2 * (5*(cos(i)^2) - 1) / (1 - e^2)^2;

            om1 = omega0 + om_dot * (t1 - t0);

            [rr, vv] = obj.convert_kep2car(a, e, i, OM1, om1, nu1, mu);
        
        
                
        end

        function [nu1, OM1, rr, vv] = propagate_orbit_increment_keplar_newton_hyperbolic(obj, a, e, i, Omega0, omega0, nu0, mu , t0, t1, R_E, J_2)

            %%%%%%%Author: Kolja Westphal, ALL RIGHTS RESERVED%%%%%%%%%%%%
            
            %%%% Input 
            
            % a [1x1] semi-major axis [km]
            % e [1x1] eccentricity [-]
            % i [1x1] inclination [rad]
            % Omega0 [1x1] RAAN [rad]
            % omega0 [1x1] argument of periapsis [rad]
            % nu0 [1x1] true anomaly [rad]
            % mu [1x1] gravitational parameter [km^3/s^2]
            % t0 [1x1] initial time [s]
            % t1 [1x1] time of progation [s]
            % R_E [1x1] body's equatorial radius [km]
            % J_2 [1x1] body's second dynamic form factor [-]
            
            %%%% Output
            
            % nu1 [1x1] propagated true anomaly
            % OM1 [1x1] propagated RAAN with J_2
            % rr [3x1] propagated position vector
            % vv [3x1] propagated velocity vector
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
            % Calculate initial values
        
            H0 = asinh( sqrt(e^2 - 1) * sin(nu0) / (1 + e*cos(nu0)) );
            M0 = e * sinh(H0) - H0;
            n = sqrt(mu/(-a^3));
            M = M0 + n * (t1 - t0);
        
            % Define the function for Newton's method
            f = @(H) M - e*sinh(H) + H;
            df = @(H) -1 * e*cosh(H) + 1;       % Changed here to use normal newton rapson method
        
            % Initial guess for E
            if e < 1.6
                if (-pi < M && M < 0) || (M > pi)
                    H_init = M - e;
                else
                    H_init = M + e;
                end
            else
                if (e < 3.6 && abs(M) > pi)
                    H_init = M - sign(M) * e;
                else
                    H_init = M / (e - 1);
                end
            end
        
            % Run Newton's method and calcualte true anomaly

            opts.return_all = true;
        
            [H, H_all] = obj.perform_newtons_method(f, df, H_init, opts);

            % True anomaly
            nu1 = atan(sqrt((e+1) / (e-1)) * tanh(H/2)) * 2;

            % For hyperbolic orbits, node drift is excluded
            OM1 = Omega0;
            om1 = omega0;

            [rr, vv] = obj.convert_kep2car(a, e, i, OM1, om1, nu1, mu);
        
                
        end
        
        function [tt, R, V, nunu, OmegaOmega] = propagate_orbit_keplar_newton(obj, a, e, i, Omega0, omega, nu0, mu , t0, t1, t_step, R_E, J_2, hyperbolic)

            %%%%%%%Author: Kolja Westphal, ALL RIGHTS RESERVED%%%%%%%%%%%%
            
            %%%% Input 
            
            % a [1x1] semi-major axis [km]
            % e [1x1] eccentricity [-]
            % i [1x1] inclination [rad]
            % Omega0 [1x1] RAAN [rad]
            % omega0 [1x1] argument of periapsis [rad]
            % nu0 [1x1] true anomaly [rad]
            % mu [1x1] gravitational parameter [km^3/s^2]
            % t0 [1x1] initial time [s]
            % t1 [1x1] time of progation [s]
            % t_step [1x1] time step of propagation [s]
            % R_E [1x1] body's equatorial radius [km]
            % J_2 [1x1] body's second dynamic form factor [-]
            % hyperbolic [bool] true for propgation of hyperbolic orbits, if not: false

            %%%% Output
            
            % tt [nx1] time array [s]
            % R [nx3] matrix of position vectors [km]
            % V [nx3] matrix of velcoty vecotrs [km]

            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                
            % Allocate memory for outputs

            Nt = floor((t1 - t0)/t_step) + 1;
            tt = linspace(t0, t1, Nt).';
        
            nunu = zeros(Nt, 1);
            OmegaOmega = zeros(Nt, 1);
            R = zeros(Nt, 3);
            V = zeros(Nt, 3);

            % Get initial position of satelite
            [rr, vv] = obj.convert_kep2car(a, e, i, Omega0, omega, nu0, mu);
            nunu(1) = nu0; % true anomaly vector
            OmegaOmega(1) = Omega0; % RAAN vector
            R(1, :) = rr'; % Position vector Matrix
            V(1, :) = vv'; % Velocity vector Matrix

            % Start propagation
            for n = 2:Nt
                t_current = tt(n);

                if hyperbolic == false
                    [nu1, Omega1, rr, vv] = obj.propagate_orbit_increment_keplar_newton(a, e, i, Omega0, omega, nu0, mu, t0, t_current, R_E, J_2);
                else
                    [nu1, Omega1, rr, vv] = obj.propagate_orbit_increment_keplar_newton_hyperbolic(a, e, i, Omega0, omega, nu0, mu, t0, t_current, R_E, J_2);
                end

        
                % Compute the orbital parameters (nu, OM, rr, vv) for the current step
                
                nunu(n) = nu1; % true anomaly vector
                OmegaOmega(n) = Omega1; % RAAN vector
                R(n, :) = rr'; % Position vector Matrix
                V(n, :) = vv'; % Velocity vector Matrix
        
                nu0 = nu1; % Update true anomaly for the next iteration
                Omega0 = Omega1;

                t0 = t_current;
            end
                
        end
    
        function [tt, R, V] = propagate_orbit_EoM_rk4(obj, rr0, vv0, t0, t1, t_step, mu, R_E, J_2)

            %%%%%%%Author: Kolja Westphal, ALL RIGHTS RESERVED%%%%%%%%%%%%
            
            %%%% Input 
            
            % rr0 [3x1] initial position vector [km]
            % vv0 [3x1] inital velocity vector [km/s]
            % t0 [1x1] initial time
            % t1 [1x1] time of progation
            % t_step [1x1] time step of propagation [s]
            % mu [1x1] gravitational parameter [km^3/s^2]
            % R_E [1x1] body's equatorial radius [km]
            % J_2 [1x1] body's second dynamic form factor [-]
            
            %%%% Output
            
            % tt [nx1] time array [s]
            % R [nx3] matrix of position vectors [km]
            % V [nx3] matrix of velcoty vecotrs [km] 
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                
            % Define derivate functions

            f = @(t, r, z) z;
            g = @(t, rr, zz) ...
                [
                -(mu*rr(1)) / (norm(rr)^3) * (1 - (3/2)*J_2 * ((R_E/norm(rr))^2 * (5*((rr(3)^2)/norm(rr)^2) - 1)));
                -(mu*rr(2)) / (norm(rr)^3) * (1 - (3/2)*J_2 * ((R_E/norm(rr))^2 * (5*((rr(3)^2)/norm(rr)^2) - 1)));
                -(mu*rr(3)) / (norm(rr)^3) * (1 - (3/2)*J_2 * ((R_E/norm(rr))^2 * (5*((rr(3)^2)/norm(rr)^2) - 3)));
                ];

            [tt, R, V] = obj.perform_rk4_2nd_ODE(f, g, t0, rr0, vv0, t1, t_step);

        end
                
        function [tt, R, V] = propagate_rk4_low_thrust(...
            obj, r0, v0, t0, tf, dt, mu, amax, ck, deltak)

            %%%%%%%Author: Tristan De La Cruz Hachiles, TUB 2026, ALL RIGHTS RESERVED%%%%%%%%%%%%
            
            %%%%%% Inputs:

            % obj - Object instance (contains auxiliary methods like convert_car2kep)
            % r0 - Initial position vector [3x1] in km
            % v0 - Initial velocity vector [3x1] in km/s
            % t0 - Initial time [s]
            % tf - Final time [s]
            % dt - Time step for RK4 integration [s]
            % mu - Gravitational parameter [km^3/s^2]
            % amax - Maximum thrust acceleration magnitude [km/s^2]
            % ck - Throttle coefficient (0 = no thrust, 1 = full thrust)
            % deltak - Thrust direction angle relative to velocity plane [rad]
            
            %%%% Outputs:

            % tt - Time history [Nx1] in seconds
            % R - Position history [Nx3] in km
            % V - Velocity history [Nx3] in km/s

            N = ceil((tf - t0)/dt);

            tt = zeros(N+1,1);
            R  = zeros(N+1,3);
            V  = zeros(N+1,3);

            r = r0(:);
            v = v0(:);

            tt(1) = t0;
            R(1,:) = r';
            V(1,:) = v';

            for i = 1:N
                % --- RK4 ---
                [k1r, k1v] = dynamics_(r, v);
                [k2r, k2v] = dynamics_(r + 0.5*dt*k1r, v + 0.5*dt*k1v);
                [k3r, k3v] = dynamics_(r + 0.5*dt*k2r, v + 0.5*dt*k2v);
                [k4r, k4v] = dynamics_(r + dt*k3r,     v + dt*k3v);

                r = r + dt/6 * (k1r + 2*k2r + 2*k3r + k4r);
                v = v + dt/6 * (k1v + 2*k2v + 2*k3v + k4v);

                tt(i+1) = tt(i) + dt;
                R(i+1,:) = r';
                V(i+1,:) = v';
            end

            % ===============================
            function [rdot, vdot] = dynamics_(r, v)
                % Computes derivatives for RK4

                %%%%%%%Author: Tristan De La Cruz Hachiles, TUB 2026, ALL RIGHTS RESERVED%%%%%%%%%%%%
                %%%% Inputs:

                % r - position vector [3x1]
                % v - velocity vector [3x1]
                
                %%%% Outputs:

                % rdot - time derivative of position = velocity [3x1]
                % vdot - time derivative of velocity = acceleration [3x1]

                rnorm = norm(r);

                % Gravity (two-body)
                a_grav = -mu * r / rnorm^3;

                % Thrust acceleration
                if ck > 0
                    t_hat = v / norm(v);
                    h = cross(r,v);
                    n_hat = h / norm(h);

                    a_thrust = ck * amax * ...
                        (cos(deltak)*t_hat + sin(deltak)*n_hat);
                else
                    a_thrust = [0;0;0];
                end

                rdot = v;
                vdot = a_grav + a_thrust;
            end
            % ===============================

        end

        function J = objective_low_thrust(obj, z, params)
            %%%%%%%Author: Tristan De La Cruz Hachiles, TUB 2026, ALL RIGHTS RESERVED%%%%%%%%%%%%
            
            % OBJECTIVE_LOW_THRUST Computes a penalty-only objective function for low-thrust GEO rendezvous.
            %
            %%%% Inputs:

            % obj - Object instance (contains propagate_rk4_low_thrust and convert_car2kep)
            % z - Decision variable vector [2N x 1], first N = delta angles, next N = throttle coefficients
            % params - Struct with fields:
                % N - Number of thrust arcs
                % r0, v0 - Initial position and velocity [3x1]
                % tw - Duration of coast phase [s]
                % dt - Time step [s]
                % mu - Gravitational parameter [km^3/s^2]
                % amax - Maximum thrust acceleration [km/s^2]
                % tLT - Total thrust phase duration [s]
                % lambda_t0 - Initial target longitude [rad]
                % nGEO - Mean motion of GEO target [rad/s]
                % rGEO, e_GEO, i_GEO - Target orbital elements
                % wa, we, wi, wl, wu - Penalty weights
            
            %%%% Output:

            % J - Scalar penalty value representing deviation from desired GEO intercept

            N = params.N;
            delta = z(1:N);
            c     = z(N+1:2*N);

            % Initial state
            r = params.r0;
            v = params.v0;

            % 1) Coast phase
            [~, R1, V1] = obj.propagate_rk4_low_thrust( ...
                r, v, 0, params.tw, params.dt, ...
                params.mu, params.amax, 0, 0);

            r = R1(end,:)';
            v = V1(end,:)';

            % 2) Thrust arcs
            dt_arc = params.tLT / N;

            for k = 1:N
                [~, Rk, Vk] = obj.propagate_rk4_low_thrust( ...
                    r, v, 0, dt_arc, params.dt, ...
                    params.mu, params.amax, c(k), delta(k));

                r = Rk(end,:)';
                v = Vk(end,:)';
            end

            tf = params.tw + params.tLT;

            % 3) Orbital elements at tf
            [a, e, inc, ~, ~] = obj.convert_car2kep(r, v, params.mu);

            % 4) Longitude error
            lambda_s = atan2(r(2), r(1));
            lambda_t = params.lambda_t0 + params.nGEO * tf;
            dLambda = wrapTo2Pi(lambda_s - lambda_t);

            % 5) Penalty-only objective
            J = params.wa * ((a - params.rGEO)/1)^2 ...
            + params.we * (e / 1e-3)^2 ...
            + params.wi * (rad2deg(inc) / 0.1)^2 ...
            + params.wl * (rad2deg(dLambda) / 0.1)^2 ...
            + params.wu * sum(c.^2)/N;

            end 

        function [R_hist, V_hist, T_hist] = propagate_low_thrust_history(obj, z, params)

            %%%%%%%Author: Tristan De La Cruz Hachiles, TUB 2026, ALL RIGHTS RESERVED%%%%%%%%%%%%
            % PROPAGATE_LOW_THRUST_HISTORY Generates full trajectory history for plotting/analysis.
            %
            %%%% Inputs:

            % obj - Object instance
            % z - Decision variable vector [2N x 1], first N = delta angles, next N = throttle coefficients
            % params - Struct with fields as in objective_low_thrust
            %

            %%%% Outputs:

            % R_hist - Position history [Nx3] km
            % V_hist - Velocity history [Nx3] km/s
            % T_hist - Time history [Nx1] s

            N = params.N;

            % ---- Control extraction ----
            delta = z(1:N);
            c     = z(N+1:2*N);

            % ---- Initial state ----
            r = params.r0;
            v = params.v0;

            R_hist = [];
            V_hist = [];
            T_hist = [];

            % ---- Coast phase ----
            [tv, R, V] = obj.propagate_rk4_low_thrust( ...
                r, v, 0, params.tw, params.dt, ...
                params.mu, params.amax, 0, 0);

            R_hist = [R_hist; R];
            V_hist = [V_hist; V];
            T_hist = [T_hist; tv];

            r = R(end,:)';
            v = V(end,:)';
            t = T_hist(end);

            % ---- Thrust phase ----
            dt_arc = params.tLT / N;

            for k = 1:N
                [tv, R, V] = obj.propagate_rk4_low_thrust( ...
                    r, v, 0, dt_arc, params.dt, ...
                    params.mu, params.amax, c(k), delta(k));

                R_hist = [R_hist; R];
                V_hist = [V_hist; V];
                T_hist = [T_hist; t + tv];

                r = R(end,:)';
                v = V(end,:)';
                t = T_hist(end);
            end
        end

        function plot_low_thrust_history(obj, R_thrust, T_thrust, params_opt)
            
            %%%%%%%Author: Tristan De La Cruz Hachiles, TUB 2026, ALL RIGHTS RESERVED%%%%%%%%%%%%
            % PLOT_LOW_THRUST_HISTORY Plots spacecraft trajectory and GEO target.
            %
            %%%% Inputs:

            % obj - Object instance
            % R_thrust - Position history of chaser [Nx3] km
            % T_thrust - Time history corresponding to R_thrust [Nx1] s
            % params_opt - Struct with target GEO elements and parameters:
            % rGEO, e_GEO, i_GEO, OM_GEO, om_GEO, lambda_t0, nGEO, mu


            Nt = length(T_thrust);

            % ---- Target propagation (GEO) ----
            R_target = zeros(Nt,3);

            for k = 1:Nt
                nu_t = params_opt.lambda_t0 + params_opt.nGEO * T_thrust(k);

                [r_t, ~] = OrbitPropagation.convert_kep2car( ...
                    params_opt.rGEO, params_opt.e_GEO, ...
                    params_opt.i_GEO, params_opt.OM_GEO, ...
                    params_opt.om_GEO, nu_t, params_opt.mu);

                R_target(k,:) = r_t';
            end

            % ---- Plot ----
            figure; hold on;

            plot3(R_thrust(:,1), R_thrust(:,2), R_thrust(:,3), ...
                'b', 'LineWidth', 1.6);

            plot3(R_target(:,1), R_target(:,2), R_target(:,3), ...
                'r--', 'LineWidth', 1.6);

            % ---- Earth ----
            [xe,ye,ze] = sphere(40);
            surf(6378*xe, 6378*ye, 6378*ze, ...
                'FaceColor',[0.6 0.8 1], ...
                'EdgeColor','none', ...
                'FaceAlpha',0.3);

            axis equal; grid on;
            xlabel('x [km]');
            ylabel('y [km]');
            zlabel('z [km]');
            legend('Chaser','Target (GEO)','Earth');
            title('Low-Thrust GEO Rendezvous');
            view(35,25);
        end
        
    end

end
