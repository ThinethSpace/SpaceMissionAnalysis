classdef InterplanetaryTransfers
    methods (Static)
        function [x, x_all] = perform_secant_method(f, x0, x1, max_iterations, tolerance)
            %%%%%%%Author: Kolja Westphal, ALL RIGHTS RESERVED%%%%%%%%%%%%
            
            %%%% Input 
            
            % f [func] function to find roots of []
            % x0, x1 [1x1] initial guesses of  t []
            % max_iterations [1x1] maximum iterations
            % tolerance [1x1] stop criteria

            %%%% Output
            
            % x [1x1] root of f []

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % ----------------
            % Secant method.
            % ----------------

            % returns initial guess if it is a root of f(x)
            if f(x0) == 0
                x = x0;
                return
            end

            % Save iteratinos values for debugging
            x_all = zeros(1,max_iterations+1);

            % iteration
            % Initialize the current guess for the iteration
            x_current = x0;

            for k = 1:max_iterations
                
                % stores results in arrays
                x_all(k) = x_current;

                % If new guesses cause devition by zero, value is found
                if (f(x0) - f(x1)) == 0
                    x = x_current;
                    return 
                end
                
                % Check if tolerance is achieved
                if abs(f(x_current)) < tolerance
                    x = x_current;
                    return
                end

                % Find new x value
                x_current = x0 - f(x0) * (x0 - x1) / (f(x0) - f(x1));

                % Update current guess for the next iteration
                x0 = x1; x1 = x_current;
            end

            % If maximum iterations reached without convergence
            warning('Maximum iterations reached.');
            x = x_current;

        end
    end
    methods
        function [vv_1, vv_2, a] = solve_lamberts_problem_secant(obj, rr_1, rr_2, delta_theta, dt, mu, factors, max_iterations, tolerance)
            
            %%%%%%%Author: Kolja Westphal, ALL RIGHTS RESERVED%%%%%%%%%%%%
            
            %%%% Input 
            
            % rr_1, rr_1 [3x1] position vectors [km]
            % dt [1x1] delta t [sec]
            % mu [1x1] gravitational parameter

            %%%% Output
            
            % vv_1, vv_1 [3x1] velocity vectors [km/s]
            % a [1x1] semi major axis

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            arguments
               obj, rr_1, rr_2, delta_theta, dt, mu, factors, max_iterations, tolerance;
            end

            % Relevant norms
            n_rr_1 = norm(rr_1); n_rr_2 = norm(rr_2);
            
            % Calculating chord length and semi parameter
            c = sqrt(n_rr_1^2 + n_rr_2^2 - 2*n_rr_1*n_rr_2*cos(delta_theta));
            s = (n_rr_1 + n_rr_2 + c) / 2;

            % Get alpha and beta
            function [alpha_0, beta_0] = get_alpha_0_beta_0(a)
                alpha_0 = asin(sqrt(s/(2*a))) * 2;
                beta_0 = asin(sqrt((s-c)/(2*a))) * 2;
            end
        
            % Lambert's Equation
            function y = f(a, alpha, beta)
                y = a^(3/2) * (alpha - beta - (sin(alpha) - sin(beta))) - sqrt(mu) * dt;
            end

            % Minimum energy solution
            a_m = s / 2;   
            [alpha_0, beta_0] = get_alpha_0_beta_0(a_m);
            t_m = f(a_m, alpha_0, beta_0);

            % Determine alpha and beta based on transfer differentiation
            if (delta_theta) <= pi && (dt <= t_m)
                beta = beta_0; alpha = alpha_0;
            elseif (delta_theta) <= pi && (dt > t_m)
                beta = beta_0; alpha = 2*pi - alpha_0;
            elseif (delta_theta) > pi && (dt <= t_m)
                beta = -1 * beta_0; alpha = alpha_0;
            elseif (delta_theta) > pi && (dt > t_m)
                beta = -1 * beta_0; alpha = 2*pi - alpha_0;
            else
                error("Angle between vectors and min energy transfer do not correlate. Abort")
            end

            % Set initial guesses for a
            a_n1 = factors(1) * a_m;
            a_n2 = factors(2) * a_m;

            % Pass function
            f_pass = @(x) f(x, alpha, beta);
            
            % Perform secant method to find a
            [a, a_all] = obj.perform_secant_method(f_pass, a_n1, a_n2, max_iterations, tolerance);

            % Determine velocities
            u_1 = rr_1 / n_rr_1;
            u_2 = rr_2 / n_rr_2;
            u_c = (rr_2 - rr_1) / c;

            A = sqrt(mu/(4*a)) * cot(alpha/2);
            B = sqrt(mu/(4*a)) * cot(beta/2);

            vv_1 = (B+A) * u_c + (B-A) * u_1;
            vv_2 = (B+A) * u_c + (B-A) * u_2;


            



        end
    end

    methods

    end
end