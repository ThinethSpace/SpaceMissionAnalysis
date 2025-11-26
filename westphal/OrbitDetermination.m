classdef OrbitDetermination
    %FUNCTIONS Summary of this class goes here
    %   Detailed explanation goes here
    methods (Static)
        function [vv_1, vv_2] = solve_lamperts_problem_gauss(rr_1, rr_2, t_1,t_2, mu, max_iterations, tolerance)
            %%%%%%%Author: Kolja Westphal, TUB 2025, ALL RIGHTS RESERVED%%%%%%%%%%%%
            
            %%%% Input 
            
            % rr_1 [3x1] position vector 1 [km]
            % rr_2 [3x1] position vector 2 [km]
            % t_1 [1x1] time of rr_1 [s]
            % t_2 [1x1] t time of rr_2 [s]

            % mu [1x1] gravitational parameter [km^3/s^2]

            % Optional:
            % max_iterations [1x1] maximum iterations to solve x_1
            % tolerance [1x] stop criteria for iterations: abs(y1-y0) < tolerance
            
            %%%% Output
            
            % vv_1 [3x1] velocity vector at rr_1 [km/s]
            % vv_2 [3x1] velocity vector at rr_2 [km/s]
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            arguments
                rr_1, rr_2, t_1, t_2, mu;
                max_iterations = 200;
                tolerance = 1E-12;
            end
            % Init arrays

            vv_1 = zeros(3);
            vv_2 = zeros(3);

            % Calculate norm of vectors 
            n_rr_1 = norm(rr_1);
            n_rr_2 = norm(rr_2);

            % check if angle between the two vectors exceeds 90°
            if (dot(rr_1, rr_2) < 0)
                error('Angle between vectors exceeds 90°')
            end
            
            % Define l and m
            delta_Omega = acos(dot(rr_1, rr_2) / ((n_rr_1 * n_rr_2)));
            l = (n_rr_1 + n_rr_2) / (4*sqrt(n_rr_1*n_rr_2) * cos(delta_Omega/2)) - 1/2;
            m = (mu*(t_2 - t_1)^2) / (2*sqrt(n_rr_1*n_rr_2) * cos(delta_Omega/2))^3;

            % Loop to determine y
            % Initial guess for y
            y0 = 1;

            for k = 1:max_iterations
                x_1 = m / y0^2;
                
                % Determine x_2
                term = 1;
                x_2 = term;
                nTerms = 4; %number of terms in the series
                for i=2:nTerms
                    j=2*i+1;
                    term=term*x_1*(j+1)/j;
                    x_2=x_2+term;
                end

                % Solve for y
                y1 = 1 + x_2*(l + x_1);
                    
                % Break 
                if (abs(y1-y0) < tolerance)
                    break
                else
                    y0 = y1;
                end
            end

            % solve f, g and g'
            p = (n_rr_1 * n_rr_2 * (1 - cos(delta_Omega))) / ...
                (n_rr_1 + n_rr_2 - 2*sqrt(n_rr_1*n_rr_2) * ...
                cos(delta_Omega/2) * (1 - 2*x_1));

            f = 1 - (n_rr_2 / p) * (1 - cos(delta_Omega));
            g = (n_rr_1*n_rr_2 * sin(delta_Omega)) / sqrt(mu * p);
            g_dash = 1 - (n_rr_1/p) * (1 - cos(delta_Omega));

            % Solve v_1 and v_2

            vv_1 = (rr_2 - f*rr_1) / g;
            vv_2 = (g_dash*rr_2 - rr_1) / g;
        end
        
        function [vv_1, vv_2, vv_3] = solve_gibbs_method(rr_1, rr_2, rr_3, min_angle_coplanar, min_angle_separation)
            
            %%%%%%%Author: Kolja Westphal, TUB 2025, ALL RIGHTS RESERVED%%%%%%%%%%%%
            
            %%%% Input 
            
            % rr_1 [3x1] position vector 1 [km]
            % rr_2 [3x1] position vector 2 [km]
            % rr_2 [3x1] position vector 3 [km]
            % Optional:
            % min_angle_coplanar [1x] min angle between vectors to check for coplanarity[deg]
            % min_angle_separation [1x] min angle between vectors in plane [deg]

            
            %%%% Output
            
            % vv_1 [3x1] velocity vector at rr_1 [km/s]
            % vv_2 [3x1] velocity vector at rr_2 [km/s]
            % vv_3 [3x1] velocity vector at rr_3 [km/s]

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            arguments
                rr_1, rr_2, rr_3,
                min_angle_coplanar = 5;
                min_angle_separation = 1;

            end

            % Calculate cross products and relevant norm of vectors
            zz_12 = cross(rr_1, rr_2);
            zz_23 = cross(rr_2, rr_3);
            zz_31 = cross(rr_3, rr_1);
            n_rr_1 = norm(rr_1); n_rr_2 = norm(rr_2); n_rr_3 = norm(rr_3);

            % Check coplanarity
            alpha_cop = asin( dot(zz_23, rr_1) / (n_z23 * n_rr_1) );

            if rad2deg(alpha_cop) > min_angle_coplanar
                error("Vectors coplanarity angle exceeds " + min_angle_separation + " degrees");
            end
            
            % Check angle difference between vectors
            a12 = acos(dot(rr_1, rr_2) / (n_rr_1 * n_rr_2));
            a23 = acos(dot(rr_2, rr_3) / (n_rr_2 * n_rr_3));
            
            if (rad2deg(a12) < min_angle_separation) || (rad2deg(a23) < min_angle_separation)
                error("Angle between vectors is lower than " + min_angle_separation + " degrees");
            end

            % Calculate vectors n, d and s
            nn = n_rr_1*zz_23 + rr_2*zz_31 + rr_3*zz_12;
            dd = zz_23 + zz_31 + zz_12;
            ss = (rr_2 - rr_3)*rr_1 + (rr_3 - rr_1)*rr_2 + (rr_1 - rr_2)*rr_3;

            B = cross([dd,dd,dd], [rr_1,rr_2,rr_3]);

            L_g = sqrt(mu/(norm(nn) * norm(dd)));
            V = B .* [L_g/n_rr_1, L_g/n_rr_2, L_g/n_rr_3] + L_g*[ss,ss,ss];

            vv_1 = V(:,1);
            vv_2 = V(:,2);
            vv_3 = V(:,3);

        end
        function [vv_2] = solve_gibbs_herrik_method(rr_1, rr_2, rr_3,tt, min_angle_coplanar, min_angle_separation)
            
            %%%%%%%Author: Kolja Westphal, TUB 2025, ALL RIGHTS RESERVED%%%%%%%%%%%%
            
            %%%% Input 
            
            % rr_1 [3x1] position vector 1 [km]
            % rr_2 [3x1] position vector 2 [km]
            % rr_2 [3x1] position vector 3 [km]
            % tt [3x1] times of the three position vectors [s]
            % Optional:
            % min_angle_coplanar [1x] min angle between vectors to check for coplanarity[deg]
            % min_angle_separation [1x] min angle between vectors in plane [deg]

            
            %%%% Output
            
            % vv_2 [3x1] velocity vector at rr_2 [km/s]

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            arguments
                rr_1, rr_2, rr_3,tt,
                min_angle_coplanar = 5;
                min_angle_separation = 5;
            end

            % Calculate cross products and relevant norm of vectors
            zz_12 = cross(rr_1, rr_2);
            zz_23 = cross(rr_2, rr_3);
            zz_31 = cross(rr_3, rr_1);
            n_rr_1 = norm(rr_1); n_rr_2 = norm(rr_2); n_rr_3 = norm(rr_3);

            % Check coplanarity
            alpha_cop = asin( dot(zz_23, rr_1) / (n_z23 * n_rr_1) );

            if rad2deg(alpha_cop) > min_angle_coplanar
                error("Vectors coplanarity angle exceeds " + min_angle_separation + " degrees");
            end
            
            % Check angle difference between vectors
            a12 = acos(dot(rr_1, rr_2) / (n_rr_1 * n_rr_2));
            a23 = acos(dot(rr_2, rr_3) / (n_rr_2 * n_rr_3));
            
            if (rad2deg(a12) > min_angle_separation) || (rad2deg(a23) > min_angle_separation)
                error("Angle between vectors is higher than " + min_angle_separation + " degrees");
            end

            dt_21 = tt(2) - tt(1);
            dt_32 = tt(3) - tt(2);
            dt_31 = tt(3) - tt(1);
            
            n_rr_1 = norm(rr_1);
            n_rr_2 = norm(rr_2);
            n_rr_3 = norm(rr_3);
            
            c1 = - (dt_32) * (1 / (dt_21 * dt_31) + mu / (12 * n_rr_1^3));
            c2 =   (dt_32 - dt_21) * (1 / (dt_21 * dt_32) + mu / (12 * n_rr_2^3));
            c3 =     (dt_21) * (1 / (dt_32 * dt_31) + mu / (12 * n_rr_3^3));
            
            vv_2 = c1*rr_1 + c2*rr_2 + c3*rr_3;

        end

        function R = solve_angles_only_approach(LOS, R_GS, tt, mu)
            %%%%%%%Author: Kolja Westphal, TUB 2025, ALL RIGHTS RESERVED%%%%%%%%%%%%
            
            %%%% Input 
            
            % LOS [3x3] Line-of-sight matrix with vectors as columns [-]
            % R_GS [3x3] groundstation position matrix with vectors as columns [-]
            % tt [3x1] time vectors [s]
            % mu [1x] gravity constant []
            
            %%%% Output
            
            % R [3x3] Position vectors for LOS vectros[km]

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Calculate intermediate matrix
            M = inv(LOS)\R_GS;

            % Form coefficients
            tau_1 = tt(1) - tt(2);
            tau_3 = tt(3) - tt(2);
            a_1 = tau_3 / (tau_3 - tau_1);
            a_1u = (tau_3*((tau_3 - tau_1)^2) - tau_3^2) / (6*(tau_3 - tau_1));
            a_3 = -1 * tau_1 / (tau_3 - tau_1);
            a_3u = -1 * (tau_1*(tau_3 * tau_1)^2 - tau_1^2) / (6*(tau_2 - tau_1));

            % Calculate d_1 and d_2 and C
            d_1 = M(2,1)*a_1 - M(2,2) + M(2,3)*a_3;
            d_2 = M(2,1)*a_1u  + M(2,3)*a_3u;
            C = dot(LOS(:,2), R_GS(:,2));

            % Get roots of polynomial
            poly = [
                1, 0, -(d_1^2 + C*d_1 + norm(R_GS(:,2))^2), 0, ...
                0, 2*mu*(C*d_2 + d_2*d_1), 0, 0, mu^2*d_2^2
                ];

            rts = roots(poly);

            % Check for complex roots
            if any(imag(r) ~= 0)
                error("Polynomial has complex roots");
            end
    
            % Calculate coefficients
            u = mu/(max(rts)^3);
            c_1 = a_1 + a_1u*u;
            c_2 = -1;
            c_3 = a_3 + a_3u*u;

            % Determine slat ranges
            cc = [c_1;c_2;c_3];
            rhorho = dot(M, -1*cc) ./ cc;

            % Determine the position vectors
            R = LOS .* rhorho' + R_GS;

        end

    end
    methods
        
        function [R, vv_2] = solve_anlges_only_approach_extended(obj, LOS, R_GS, tt, mu, tol, max_iterations)
            %%%%%%%Author: Kolja Westphal, TUB 2025, ALL RIGHTS RESERVED%%%%%%%%%%%%
            
            %%%% Input 
            
            % LOS [3x3] Line-of-sight matrix with vectors as columns [-]
            % R_GS [3x3] groundstation position matrix with vectors as columns [-]
            % tt [3x1] time vectors [s]
            % mu [1x1] gravity constant []
            
            % Optional:
            % tol [1x1] Tolerance of change in slant ranges for break of iteration [-]
            % max_iterations [1x] maximum iterations of loop to determine slant ranges
            
            %%%% Output
            
            % R [3x3] Position vectors for LOS vectros[km]

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Solve angles only approach for first guess of roh -> R
            R = obj.solve_angles_only_approach(LOS, R_GS, tt, mu);

            % Init empty slant range vector
            rhorho = zeros(3,1);

            % Start the loop
            n = 1;
            while n <= max_iteratinos

                % Check angles between vectors
                a12 = acos(dot(R(:,1), R(:,2)) / (norm(R(:,1)) * norm(R(:,2))));
                a23 = acos(dot(R(:,2), R(:,3)) / (norm(R(:,2)) * norm(R(:,3))));
                
                if (rad2deg(a12) > 3) && (rad2deg(a23) > 3)
                   [~, vv_2, ~] = obj.solve_gibbs_method(R(:,1), R(:,2), R(:,3));
                elseif (rad2deg(a12) < 3) && (rad2deg(a23) < 3)
                    vv_2 = obj.solve_gibbs_herrik_method(R(:,1), R(:,2), R(:,3));
                else
                    error('No method could be applied to solve for v_2')
                end
                
                u_dot = -1*mu*(R(:,2) * vv_2 / norm(R(:,2))) / norm(R(:,2))^4;
                tau_1 = tt(1) - tt(2);
                tau_3 = tt(3) - tt(2);
                u = mu/(norm(R(:,2))^3);

                f = @(tau_i) 1 - (1/2)*u*tau_i^2 - (1/6)*u_dot*tau_i^3 - ...
                    (1/24)*(u^2)*tau_i^4 - (1/30)*u*u_dot*tau_i^5;
                g = @(tau_i) tau_i - (u/6)*tau_i^3 - (u_dot/12)*tau_i^4 - ...
                    ((u^2)/120)*tau_i^5 - (u*u_dot/120)*tau_i^6;

                c_1 = g(tau_3) / (f(tau_1)*g(tau_3) - f(tau_3)*g(tau_1));
                c_2 = -1;
                c_3 = +1*g(tau_1) / (f(tau_1)*g(tau_3) - f(tau_3)*g(tau_1));

                cc = [c_1;c_2;c_3];
                rhorho_new = dot(M, -1*cc) ./ cc;

                R = LOS .* rohroh' + R_GS;

                % Check if the change in the slat ranges is almost zero
                if norm(rhorho - rhorho_new) < tol
                    break;
                else
                    rhorho = rhorho_new;
                    n = n + 1;
                end

            end
        end
            
    end

end
