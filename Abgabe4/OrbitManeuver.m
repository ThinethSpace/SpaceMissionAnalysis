classdef OrbitManeuver
    methods (Static)
        function [v_initial, v_trans_a, v_trans_b, v_final, phi_trans_b, t_trans] = one_tangent_burn_LEO_GEO(r_initial, r_final, nu_trans_b, mu)
            
                R_invert = r_initial / r_final;
                e_trans = (R_invert - 1) / (cos(nu_trans_b) - R_invert);
                a_trans = r_initial / (1 - e_trans);

                if (r_initial/2 > a_trans) || (r_final/2 > a_trans)
                    t_trans = 0;
                    v_initial = 0; v_trans_a = 0; v_trans_b = 0; v_final = 0; phi_trans_b = 0;
                    return;
                end                    

                v_initial = sqrt(mu / r_initial);
                v_final = sqrt(mu / r_final);
                v_trans_a = sqrt((2*mu / r_initial) - (mu/a_trans));
                v_trans_b = sqrt((2*mu / r_final) - (mu/a_trans));

                %delta_v_a = v_trans_a - v_initial;
                phi_trans_b = atan((e_trans*sin(nu_trans_b)) / (1 + e_trans * cos(nu_trans_b)));
            
                %delta_v_b = sqrt(v_trans_b^2 + v_final^2 - 2*v_trans_b*v_final*cos(phi_trans_b));

                E = acos((e_trans + cos(nu_trans_b)) / (1 + e_trans*cos(nu_trans_b)));

                % Simplified version of time of flight eqt because E0 = 0
                t_trans = sqrt(a_trans^3 / mu) * (E - e_trans*sin(E));
            end
    end
    methods
        function [delta_v_phase, delta_v_trans1, delta_v_trans2, t_total] = perform_noncoplanar_phasing_circular_orbits( ...
            obj, a_initial, i_initial, RAAN_initial, a_target, i_target, lambda_target_t0, u, mu, k_tgt, split)
            %%%%%%%Author: Kolja Westphal, ALL RIGHTS RESERVED%%%%%%%%%%%%
            
            %%%% Input 

            %%%% Output
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Determine relevant parameVters for transfer orbit (Hohmann Transfer)
            % Even though, the transfer doesn't happen from the initial orbit,
            % because the transfer starts at perigee, r_initial is the same for initial and phase orbit
            r_initial = a_initial;
            r_target = a_target;
            a_transfer = (r_initial + r_target) / 2;
            t_transfer = pi * sqrt(a_transfer^3 / mu);

            % Calculate angular velocity for initial and target orbit
            omega_initial = sqrt(mu / r_initial^3);
            omega_target = sqrt(mu / r_target^3);

            % Lead angle at t2
            alpha_lead_t2 = omega_target * t_transfer;

            % Angle from initial position of interceptor to node
            % TODO: this equation can also be 2*pi - u, depending if the transfer happens at the ascending or decending node
            delta_phase_angle_initial = pi - u;

            % Time from initial point t0 of interceptor on initial orbit to node
            delta_t_node = delta_phase_angle_initial / omega_initial;

            % True longitude at t1 of target when interceptor is at node
            lambda_true_t1 = lambda_target_t0 + omega_target * delta_t_node;

            % Phase angle at t1 between target and interceptor 
            phase_angle_t1 = (RAAN_initial + pi) - lambda_true_t1;

            % Lead angle at t1
            alpha_lead_t1 = pi + phase_angle_t1;

            % Determine phase orbit to time the transfer
            rho_phase = ((alpha_lead_t1 - alpha_lead_t2) + 2*pi*k_tgt) / omega_target;
            k_int = k_tgt + 1;
            a_phase = (mu*(rho_phase / (k_int * 2*pi))^2)^(1/3);

            % Relevant velocities
            v_initial = sqrt(mu/a_initial);
            v_target = sqrt(mu/a_target);
            v_phase = sqrt(2*mu / a_initial - mu/a_phase);
            v_trans1 = sqrt(2*mu / a_initial - mu/a_transfer);
            v_trans2 = sqrt(2*mu / a_target - mu/a_transfer);

            t_total = delta_t_node + 2*pi*sqrt((a_phase^3) / mu) + t_transfer;

            if ~split
                % Delta vs
                delta_v_phase = abs(v_phase - v_initial);
                delta_v_trans1 = abs(v_trans1 - v_phase);
                delta_v_trans2 = sqrt(v_trans2^2 + v_target^2 - 2*v_trans2*v_target*cos(i_initial - i_target));
            else

                % calculate s-ratio with Lisowski method
                factor1 = 1 / (i_initial - i_target);
                factor2 = (v_initial * v_trans1) / (v_target * v_trans2);

                % Liswowski s-ratio
                delta_i = i_initial-i_target;
                s = factor1 * atan( sin(delta_i) ... 
                    / (factor2 + cos(delta_i)) ); % rad
                % s = 0.006319440565421 % ARIS' numerical result

                % total deltaV for two burns at GEO
                delta_v_phase = abs(v_phase - v_initial);
                delta_v_trans1 = sqrt(v_initial ^ 2 + v_trans1 ^ 2 - 2 * v_initial * v_trans1 * cos(s * delta_i)); % km/s
                delta_v_trans2 = sqrt(v_target ^ 2 + v_trans2 ^ 2 - 2 * v_target * v_trans2 * cos((1 -s) * delta_i)); % km/s
            end


            
        end
        function [x_opt_all, t_w_opt_all, eta_opt_all] = perform_three_impulse_transfer(...
            obj, a_initial, i_initial, RAAN_initial, a_target, i_target, lambda_target_t0, u_t0, mu, allowed_angle, num_node_crossings)
            %%%%%%%Author: Kolja Westphal, ALL RIGHTS RESERVED%%%%%%%%%%%%
            
            %%%% Input 

            %%%% Output
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % For the this function, we assume a 3 burn transfer
            % The first burn is an combined burn
            % The second burn is a circulization burn to reach geo height
            % The last burn is a pure inclination burn

            % Determine relevant parameVters for transfer orbit (Hohmann Transfer)
            r_initial = a_initial;
            r_target = a_target;
            a_transfer = (r_initial + r_target) / 2;
            t_transfer = pi * sqrt(a_transfer^3 / mu);

            % Calculate angular velocity for initial and target orbit
            omega_initial = sqrt(mu / r_initial^3);
            omega_target = sqrt(mu / r_target^3);

            function delta_lambda = determine_delta_lambda(t_w)

                % Seperation of interceptor and target after 3rd burn
                
                lambda_s1 = mod(RAAN_initial + u_t0 + omega_initial * t_w, 2*pi);
                angle = pi - mod(RAAN_initial + u_t0 + omega_initial * t_w, pi);
                lambda_s2 = mod(lambda_s1 + pi, 2*pi);
                lambda_s3 = mod(lambda_s2 + angle, 2*pi);
                
                lambda_t1 = mod(lambda_target_t0 + omega_target * t_w, 2*pi);
                lambda_t2 = mod(lambda_t1 + omega_target * t_transfer, 2*pi);
                lambda_t3 = mod(lambda_t2 + angle, 2*pi);

                %x_s = r_target*cos(lambda_s);
                %x_tgt = r_target*cos(lambda_t);

                %y_s = r_target*sin(lambda_s);
                %y_tgt = r_target*sin(lambda_t);

                %lambda_t_final = atan2(y_tgt, x_tgt);
                %lambda_s_final = atan2(y_s, x_s);

                delta_lambda = abs(wrapToPi(lambda_s3 - lambda_t3));

            end

            x = linspace(1, 3.5E4, 10000);
            y = arrayfun(@determine_delta_lambda, x);
            plot(x, y)

            function delta_v = determine_delta_v(eta)
                
                % Relevant velocities
                v_initial = sqrt(mu/a_initial);
                v_target = sqrt(mu/a_target);
                v_trans1 = sqrt(2*mu / a_initial - mu/a_transfer);

                delta_v1 = sqrt(v_initial^2 + v_trans1^2 - 2*v_initial*v_trans1*cos(eta * (i_initial - i_target)));
                delta_v2 = v_target - sqrt(mu*((2/a_target) - (1/a_transfer)));
                delta_v3 = 2*sqrt(mu/a_target) * sin(((1-eta)*(i_initial - i_target)) / 2);

                delta_v = delta_v1 + delta_v2 + delta_v3;


            end

            % Determine box constrains
            % Determine time of first node crossing (ascending or decending node)
            t_cross1 = omega_initial * min(pi - u_t0, 2*pi - u_t0);

            % Allowed times of transfer close to a node
            tt_cross1_box = [t_cross1 - omega_initial * allowed_angle, t_cross1 + omega_initial * allowed_angle];

            % Determine period on inital orbit
            period = 2*pi*sqrt((a_initial^3) / mu);

            box_constraints = zeros(num_node_crossings, 2);
            for i = 1:num_node_crossings
                box_constraints(i,:) = tt_cross1_box + i*period/2;
            end

            % Optimize in box constraints
            w_lambda = 1;  % weighting factor

            % Preallocate results
            x_opt_all = zeros(num_node_crossings, 1)
            t_w_opt_all = zeros(num_node_crossings,1);
            eta_opt_all = zeros(num_node_crossings,1);
                    
            for i = 1:num_node_crossings
                % Determine time window for this node crossing
                t_lower = box_constraints(i, 1);
                t_upper = box_constraints(i, 2);
                
                % Initial guesses
                t0 = (t_lower + t_upper)/2;
                eta0 = 0.5;  % example initial guess for eta
                
                % Combine variables for optimization: x = [t_w, eta]
                x0 = [t0, eta0];
                
                % Lower and upper bounds
                lb = [t_lower, 0];    % eta >= 0
                ub = [t_upper, 1];    % eta <= 1
                
                % Define objective function
                J_fun = @(x) determine_delta_v(x(2)) + w_lambda * determine_delta_lambda(x(1))^2;
                
                % Run fmincon
                options = optimoptions('fmincon','Display','off');
                [x_opt, fval] = fmincon(J_fun, x0, [], [], [], [], lb, ub, [], options);
                
                % Store results
                x_opt_all(i) = fval;
                t_w_opt_all(i) = x_opt(1);
                eta_opt_all(i) = x_opt(2);
            end





        end
        function [x_opt_all, t_w_opt_all, delta_v_all, delta_lambda_all, nu_opt_all] = perform_two_impulse_transfer(obj, a_initial, i_initial, RAAN_initial, a_target, i_target, lambda_target_t0, u_t0, mu, num_half_orbit_revs)

            % Determine relevant parameters for transfer orbit
            r_initial = a_initial;
            r_target = a_target;

            % Calculate angular velocity for initial and target orbit
            omega_initial = sqrt(mu / r_initial^3);
            omega_target = sqrt(mu / r_target^3);

            % Angle to node for interceptor
            nu = @(t_w) pi - mod(u_t0 + omega_initial * t_w, pi);

            t_w_allowed = @(nu_val, k) (pi - nu_val + k*pi - u_t0) / omega_initial;

            %x = linspace(1, 3.5E3, 1000);
            %%y = arrayfun(@(t) one_tangent_burn(r_initial, r_target, phi(t)), x);
            %y = arrayfun(@(t) one_tangent_burn(r_initial, r_target, nu(t)), x);
            %y2 = arrayfun(@(t) nu(t), x);
            %figure
            %yyaxis left
            %plot(x, y)
            %yyaxis right
            %plot(x, y2)

            function delta_v = determine_delta_v(v_initial, v_trans_a, v_trans_b, v_final, phi_trans_b)

                % Tangential burn im Leo
                delta_v_a = v_trans_a - v_initial;

                %Circulization burn im Geo
                delta_v_b = sqrt(v_final^2 + v_trans_b^2 - 2*v_final*v_trans_b*cos(phi_trans_b));

                % Inclination change in Geo
                delta_v_inc = 2*v_final*sin((i_initial - i_target) / 2);

                delta_v = delta_v_a + delta_v_b + delta_v_inc;

            end

            function delta_lambda = determine_delta_lambda(t_w, t_transfer)                
                delta_lambda = abs(wrapToPi((RAAN_initial + pi) - lambda_target_t0 + omega_target * (t_w + t_transfer)));
            end

            % Objective function

            function val = J(t_w, w_lambda)
                nu_current = nu(t_w);
                [v_initial, v_trans_a, v_trans_b, v_final, phi_trans_b, t_trans] = obj.one_tangent_burn_LEO_GEO(r_initial, r_target, nu_current, mu);

                val = determine_delta_v(v_initial, v_trans_a, v_trans_b, v_final, phi_trans_b) + ...
                w_lambda * determine_delta_lambda(t_w, t_trans);
            end

            function [delta_v, delta_lambda, nu_current] = J_sol(t_w)
                nu_current = nu(t_w);
                [v_initial, v_trans_a, v_trans_b, v_final, phi_trans_b, t_trans] = obj.one_tangent_burn_LEO_GEO(r_initial, r_target, nu_current, mu);

                delta_v = determine_delta_v(v_initial, v_trans_a, v_trans_b, v_final, phi_trans_b);
                delta_lambda = determine_delta_lambda(t_w, t_trans);
            end

            %%%%% Optimization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Preallocate results
            x_opt_all = zeros(num_half_orbit_revs, 1);
            t_w_opt_all = zeros(num_half_orbit_revs,1);

            delta_v_all = zeros(num_half_orbit_revs,1);
            delta_lambda_all = zeros(num_half_orbit_revs,1);
            nu_opt_all = zeros(num_half_orbit_revs,1);

            % Start loop
            for i = 1:num_half_orbit_revs

                % Determine time window for this node crossing
                t_lower = t_w_allowed(pi, i);
                t_upper = t_w_allowed(2.4, i);
                
                % Initial guesses
                t0 = (t_lower + t_upper)/2;
                x0 = t0;
                
                % Lower and upper bounds
                lb = t_lower;
                ub = t_upper;

                % Weighting factor
                w_lambda = 1E1;

                % Define objective function
                J_fun = @(x) J(x, w_lambda);
                
                % Run fmincon
                options = optimoptions('fmincon','Display','off');
                [x_opt, fval] = fmincon(J_fun, x0, [], [], [], [], lb, ub, [], options);
                
                % Store results
                x_opt_all(i) = fval;
                t_w_opt_all(i) = x_opt;

                [delta_v_all(i), delta_lambda_all(i), nu_opt_all(i)] = J_sol(x_opt);
            end
            

        end
    end
end