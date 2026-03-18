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
            x = NaN;

        end
        
       function T = parse_horizon_file(filename, output_year)
            % Parse JPL Horizons output file (position + velocity)
            % Input:  filename - path to Horizons text file
            % Output: 
            % T [table] Data table with dates, position, vecocity
            fid = fopen(filename,'r');
            if fid == -1
                error('Cannot open file: %s', filename);
            end
    
            dates = datetime.empty(0,1); % initialize as datetimedates = [];
            r = [];
            v = [];
    
            in_data = false;
    
            while ~feof(fid)
                line = strtrim(fgetl(fid));

                % Start of data
                if contains(line,'$$SOE')
                    in_data = true;
                    continue
                end

                % End of data
                if contains(line,'$$EOE')
                    break
                end

                if in_data
                    tokens = strsplit(line,',');
                    if numel(tokens) < 8
                        continue
                    end

                    % Parse calendar date (TDB)
                    date_str = strtrim(tokens{2}); % e.g., 'A.D. 2030-Jun-01 00:00:00.0000'
                    date_str = strrep(date_str,'A.D. ',''); % remove prefix
                    dates(end+1,1) = datetime(date_str,'InputFormat','yyyy-MMM-dd HH:mm:ss.SSSS');

                    % Positions as row
                    r(end+1,:) = [str2double(tokens{3}), str2double(tokens{4}), str2double(tokens{5})];

                    % Velocities as row
                    v(end+1,:) = [str2double(tokens{6}), str2double(tokens{7}), str2double(tokens{8})];
                end
            end

            % Convert time
            if output_year == "year"
                v = v * 365.24217;
            end

            % Create table from data
            T = table(dates, r, v, 'VariableNames', {'Date', 'Position', 'Velocity'});
    
            fclose(fid);
        end

        function angle = angle_between_vectors(a, b, r)
            %%%%%%%Author: Kolja Westphal, ALL RIGHTS RESERVED%%%%%%%%%%%%
            
            %%%% Input 
            
            % a, b [3x1] vectors where angle shall be found []
            % r [3x1] rotation vectors which defines the rotation direction for the angle []


            %%%% Output
            
            % angle [1x1] angle between vectors [rad]

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %a = a / norm(a);
            %b = b / norm(b);
            %r = r / norm(r);

            % Signed angle
            %angle = atan2(dot(r, cross(a, b)), dot(a, b));
            cos_angle = dot(a,b) / (norm(a)*norm(b));

            angle = acos(cos_angle);


            % Optional: map to [0, 2*pi)
            %angle = mod(angle, 2*pi);

        end
    end
    methods
        function [vv_1, vv_2, a] = solve_lamberts_problem_secant(obj, rr_1, rr_2, delta_theta, dt, mu, factors, max_iterations, tolerance)
            
            %%%%%%%Author: Kolja Westphal, ALL RIGHTS RESERVED%%%%%%%%%%%%
            
            %%%% Input 
            
            % rr_1, rr_1 [3x1] position vectors [km]
            % delta_theta [1x1] angle between positiion vectors
            % dt [1x1] delta t [sec]
            % mu [1x1] gravitational parameter
            % factors [2x1] factors for a_m for initial guesses of a
            % max_iterations [1x1] maximum iterations for lambert solver
            % tolerance [1x1] tolerance value for lambert solver

            %%%% Output
            
            % vv_1, vv_1 [3x1] velocity vectors [km/s]
            % a [1x1] semi major axis

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            arguments
               obj, rr_1, rr_2, delta_theta, dt, mu, factors, max_iterations, tolerance;
            end

            % Relevant norms
            n_rr_1 = norm(rr_1); n_rr_2 = norm(rr_2);
            
            % Calculating chord length and semi parameter constants
            c = sqrt(n_rr_1^2 + n_rr_2^2 - 2*n_rr_1*n_rr_2*cos(delta_theta));
            s = (n_rr_1 + n_rr_2 + c) / 2;

            % Get alpha and beta
            function [alpha, beta] = get_alpha_beta(a)
                alpha = asin(sqrt(s/(2*a))) * 2;
                beta = asin(sqrt((s-c)/(2*a))) * 2;
            end

            function [alpha, beta] = correct_angles_for_transfer_type(alpha, beta, dt, t_m)
                % Determine alpha and beta based on transfer differentiation
                if (delta_theta) <= pi && (dt <= t_m)
                    beta = beta; alpha = alpha;
                elseif (delta_theta) <= pi && (dt > t_m)
                    beta = beta; alpha = 2*pi - alpha;
                elseif (delta_theta) > pi && (dt <= t_m)
                    beta = -1 * beta; alpha = alpha;
                elseif (delta_theta) > pi && (dt > t_m)
                    beta = -1 * beta; alpha = 2*pi - alpha;
                else
                    error("Angle between vectors and min energy transfer do not correlate. Abort")
                end
            end
        
            % Lambert's Equation
            function y = f(a, dt, t_m, mu)
                [alpha, beta] = get_alpha_beta(a);
                [alpha, beta] = correct_angles_for_transfer_type(alpha, beta, dt, t_m);

                y = a^(3/2) * (alpha - beta - (sin(alpha) - sin(beta))) - sqrt(mu) * dt;
            end

            % Minimum energy solution
            a_m = s / 2;
            [alpha_0_unc, beta_0_unc] = get_alpha_beta(a_m);

            if delta_theta <= pi
                beta_m = beta_0_unc;
            else
                beta_m = -beta_0_unc;
            end

            t_m = sqrt(a_m^3 / mu) * (pi - beta_m + sin(beta_m));

            % Set initial guesses for a
            a_n1 = factors(1) * a_m;
            a_n2 = factors(2) * a_m;

            f_pass = @(x) f(x, dt, t_m, mu);
            
            % Perform secant method to find a
            [a, a_all] = obj.perform_secant_method(f_pass, a_n1, a_n2, max_iterations, tolerance);

            if isnan(a)
                vv_1 = [];
                vv_2 = [];
                a = NaN;
                return 
            end
            % Get alphas and betas with iterated as
            [alpha, beta] = get_alpha_beta(a);
            [alpha, beta] = correct_angles_for_transfer_type(alpha, beta, dt, t_m);

            % Determine velocities
            u_1 = rr_1 / n_rr_1;
            u_2 = rr_2 / n_rr_2;
            u_c = (rr_2 - rr_1) / c;

            A = sqrt(mu/(4*a)) * cot(alpha/2);
            B = sqrt(mu/(4*a)) * cot(beta/2);

            vv_1 = (B+A) * u_c + (B-A) * u_1;
            vv_2 = (B+A) * u_c - (B-A) * u_2;


        end
        
        function create_porkchop_plot(obj, lambert_solver_parameters, departure_data, arrival_data)

            %%%%%%%Author: Kolja Westphal, ALL RIGHTS RESERVED%%%%%%%%%%%%
            
            %%%% Input 
            
            % lambert_solver_parameters [struct] parameters for the lambert solver in a struct [-]

            %%%% Output
            

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Determine number of ephemeris in table
            num_ephemeris_departure = height(departure_data);
            num_ephemeris_arrival = height(arrival_data);

            lsp = lambert_solver_parameters;

            % Determine roation vector of main body for angle determination
            gravity_body_rotation_vector = cross(departure_data.Position(1,:), departure_data.Position(2,:));

            % Allocate output
            delta_v_inf_arr = zeros(num_ephemeris_departure, num_ephemeris_arrival);
            C3 = zeros(num_ephemeris_departure, num_ephemeris_arrival);

            for i = 1:num_ephemeris_departure
                for j = 1:num_ephemeris_arrival
                    dt = years(arrival_data.Date(j) - departure_data.Date(i));

                    if dt < 0
                        % Create arbitrary high number for non realistic transfer
                       delta_v_inf(i, j) = 1E10;
                       continue;
                    end
                    %delta_theta = obj.angle_between_vectors(departure_data.Position(i,:), arrival_data.Position(j,:), gravity_body_rotation_vector);

                    angle_short = obj.angle_between_vectors(departure_data.Position(i,:), arrival_data.Position(j,:));

                    angle_long = 2*pi - angle_short;

                    [vv_1_s, vv_2_s, a_s] = obj.solve_lamberts_problem_secant(departure_data.Position(i,:), arrival_data.Position(j,:), angle_short, dt, lsp.mu, lsp.factors, lsp.max_iterations, lsp.tolerance);
                    [vv_1_l, vv_2_l, a_l] = obj.solve_lamberts_problem_secant(departure_data.Position(i,:), arrival_data.Position(j,:), angle_long, dt, lsp.mu, lsp.factors, lsp.max_iterations, lsp.tolerance);

                    if isnan(a_s) || isnan(a_l)
                        continue;
                    end
                    v_inf_dep_s = vv_1_s  - departure_data.Velocity(i,:);
                    v_inf_arr_s = vv_2_s - arrival_data.Velocity(j,:);

                    v_inf_dep_l = vv_1_l  - departure_data.Velocity(i,:);
                    v_inf_arr_l = vv_2_l - arrival_data.Velocity(j,:);

                    delta_v_inf_s = norm(v_inf_dep_s) + norm(v_inf_arr_s);
                    delta_v_inf_l = norm(v_inf_dep_l) + norm(v_inf_arr_l);
                    
                    if delta_v_inf_s < delta_v_inf_l
                        v_inf_dep = v_inf_dep_s;
                        v_inf_arr = v_inf_arr_s;
                    else
                        v_inf_dep = v_inf_dep_l;
                        v_inf_arr = v_inf_arr_l;
                    end


                    delta_v_inf_arr(i, j) = norm(v_inf_arr);
                    C3(i, j) = norm(v_inf_dep)^2;

                    % Formatted output
                    %fprintf('i = %d, j = %d, angle = %d\n', i, j, rad2deg(delta_theta))
                end
            end

            % Plot

            % Create figure
            % Create porkchop plot
            figure('Name','Porkchop Plot ΔV');

            % Convert datetime axes to datenums for contour plotting
            x_dn = datenum(departure_data.Date);
            y_dn = datenum(arrival_data.Date);

            % Filled contour plot of delta_v_inf over departure/arrival dates
            AU2km  = 149597870.7;
            year2s = 365.2422 * 86400;
            AUyr2kms = AU2km / year2s;
            dv_total = (sqrt(C3) + delta_v_inf_arr) * AUyr2kms;
            contourf(x_dn, y_dn, dv_total', 1:1:20, 'LineColor', 'none'); hold on;
            grid on;
            colorbar
            ylabel(colorbar, 'Transfer energy proxy: ||v_{\infty,dep}|| + ||v_{\infty,arr}||');
            xlabel('Departure Date')
            ylabel('Arrival Date')
            title('Earth–Mars Porkchop Plot (2030–2033)')

            % Format axes to show readable calendar dates
            ax = gca;
            ax.XAxis.Exponent = 0;
            ax.YAxis.Exponent = 0;
            datetick('x','yyyy-mmm-dd','keeplimits','keepticks');
            datetick('y','yyyy-mmm-dd','keeplimits','keepticks');


        end 
    end


end