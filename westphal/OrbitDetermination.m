classdef OrbitDetermination
    %FUNCTIONS Summary of this class goes here
    %   Detailed explanation goes here
    methods (Static)
        function [vv_1, vv_2] = solve_lamperts_problem_gauss(rr_1, rr_2, mu, max_iterations, tolerance)
            %%%%%%%Author: Kolja Westphal, TUB 2025, ALL RIGHTS RESERVED%%%%%%%%%%%%
            
            %%%% Input 
            
            % rr_1 [3x1] position vector 1 [km]
            % rr_2 [3x1] position vector 2 [km]
            % mu [1x1] gravitational parameter [km^3/s^2]

            % Optional:
            % max_iterations [1x1] maximum iterations to solve x_1
            % tolerance [1x] stop criteria for iterations: abs(y1-y0) < tolerance
            
            %%%% Output
            
            % vv_1 [3x1] velocity vector at rr_1 [km/s]
            % vv_2 [3x1] velocity vector at rr_2 [km/s]
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            arguments
                rr_1, rr_2, mu;
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
            l = (n_rr_1 + n_rr_2) / (4*sqrt(rr_1*rr_2) * cos(delta_Omega/2));
            m = (mu*(t_2 - t-1)^2) / (4*sqrt(rr_1*rr_2) * cos(delta_Omega/2))^3;

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
                    return
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
    end
    methods
    end

end
