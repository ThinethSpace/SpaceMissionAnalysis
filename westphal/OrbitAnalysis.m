classdef OrbitAnalysis
    methods (Static)
        function [distance,c] = haversine(lat1, lon1, lat2, lon2, R)
            %%%%%%%Author: Kolja Westphal, ALL RIGHTS RESERVED%%%%%%%%%%%%
            
            %%%% Input
            % lat1, lon1 [vector] point 1 latitude and longitude [rad]
            % lat2, lon2 [scalar] point 2 latitude and longitude [rad]
            % R [scalar] radius of the sphere 
            
            %%%% Output
            % distance [scalar] distance between the two points 
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Haversine formula to calculate distance between two points on a sphere
            dLat = lat2 - lat1;
            dLon = lon2 - lon1;
            a = sin(dLat/2).^2 + cos(lat1) .* cos(lat2) .* sin(dLon/2).^2;
            c = 2 * atan2(sqrt(a), sqrt(1-a));
            distance = R * c;
        end
        function swath_width = get_swath_width(altitude, R_earth, boresight_angle)
            % Calculate swath width based on altitude
            r_sat = R_earth + altitude;

            % Determine gamma and choose the greater angle
            gamma = asin(r_sat * sin(boresight_angle) / R_earth);
            if gamma < pi/2
                gamma = pi - gamma;
            end

            roh = R_earth * cos(gamma) + r_sat + cos(boresight_angle);
            
            Delta = asin(roh * sin(boresight_angle) / R_earth);
            swath_width = 2 * R_earth * Delta;
        end
        function [lla, pts] = create_grid(n)
            i = (0:n-1)' + 0.5;
            phi = acos(1 - 2*i/n);
            golden_ratio = (1 + sqrt(5)) / 2;
            theta = 2 * pi * i / golden_ratio;

            x = cos(theta) .* sin(phi);
            y = sin(theta) .* sin(phi);
            z = cos(phi);
            pts = [x, y, z];

            % LLA (Latitude, Longitude, Altitude)
            lat = asind(pts(:,3));
            lon = atan2d(pts(:,2), pts(:,1));
            lla = [lat, lon, zeros(size(lat))];


        end
        function e = get_elevation_angle(R_earth, rr, central_angle)

            %%%%%%%Author: Kolja Westphal, ALL RIGHTS RESERVED%%%%%%%%%%%%
            
            %%%% Input
            % R_earth [1x1] radius of the Earth [km]
            % rr [3x1] satellite position [km]
            % central_angle [1x1] central angle between satellite and ground point [rad]
            
            %%%% Output
            % e [1x1] elevation angle [rad]
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            r_sat = norm(rr);
            e = atan((r_sat .* cos(central_angle) - R_earth) ./ (r_sat .* sin(central_angle)));
        end
end

    methods
        function mrt = get_mean_revisit_time(obj, lat_grid, lon_grid, tt, LLA_sat, dt, half_swath)
            %%%%%%%Author: Kolja Westphal, ALL RIGHTS RESERVED%%%%%%%%%%%%
            
            %%%% Input
            % lat_grid [mxn] latitude grid points [deg]
            % lon_grid [mxn] longitude grid points [deg]
            % tt [nx1] time vector [s]
            % LLA_sat [nx3] satellite lat, lon, alt [deg, deg, km]
            % dt [scalar] time step [s]

            %%%% Output
            % mrt [scalar] mean revisit time []
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Constants
            R_earth = 6371;             % km

            % Create a Global Grid (Approx 2-degree spacing) and start NorthWest
            %grid_pts = [deg2rad(lat_grid(:)), deg2rad(lon_grid(:))];
            grid_pts = [deg2rad(lat_grid), deg2rad(lon_grid)];
            num_pts = size(grid_pts, 1);

            sat_lat = deg2rad(LLA_sat(:,1));
            sat_lon = deg2rad(LLA_sat(:,2));

            % Vectorized Access Calculation
            % Pre-allocate a logical matrix: Rows = Grid Points, Cols = Time
            last_access = false(num_pts,1);
            last_rise_time = nan(num_pts,1);
            revisit_sum = zeros(num_pts,1);
            revisit_count = zeros(num_pts,1);

            h = waitbar(0,'Processing...');

            num = length(tt);

            for i = 1:num

                waitbar(i/num,h);

                dist = obj.haversine(grid_pts(:,1), grid_pts(:,2), ...
                                     sat_lat(i), sat_lon(i), R_earth);

                current_access = dist <= half_swath;

                % Detect rises (0 -> 1)
                rises = current_access & ~last_access;

                % If this is not first rise, compute revisit
                valid = rises & ~isnan(last_rise_time);
                revisit_sum(valid) = revisit_sum(valid) + ...
                                     (tt(i) - last_rise_time(valid));
                revisit_count(valid) = revisit_count(valid) + 1;

                % Update rise times
                last_rise_time(rises) = tt(i);

                % Update state
                last_access = current_access;
            end

            close(h);

            mrt = revisit_sum ./ revisit_count;
            mrt(revisit_count == 0) = NaN;         

        end
        function get_mean_revisit_time_increment(obj, t_start, t_current, rr, lat_grid, lon_grid,half_swath, data) 
            %%%%%%%Author: Kolja Westphal, ALL RIGHTS RESERVED%%%%%%%%%%%%
            
            %%%% Input
            % t_start [1x1] start time of the simulation [datetime]
            % t_current [1x1] current time in the simulation [s]
            % rr [3x1] satellite position in ECI [km]
            % lat_grid [mxn] latitude grid points [deg]
            % lon_grid [mxn] longitude grid points [deg]
            % half_swath [scalar] half swath width [km]
            % data [SharedData] struct for storing intermediate results
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Constants
            R_earth = 6371;             % km

            % Create a Global Grid (Approx 2-degree spacing) and start NorthWest
            %grid_pts = [deg2rad(lat_grid(:)), deg2rad(lon_grid(:))];
            grid_pts = [deg2rad(lat_grid), deg2rad(lon_grid)];

            % Get current LLA
            t_current_seconds = seconds(t_current);
            time = t_start + t_current_seconds;
            sat_lla = eci2lla(rr'*1000, datevec(time));
            sat_lat = deg2rad(sat_lla(1));
            sat_lon = deg2rad(sat_lla(2));

            % Vectorized Access Calculation
            dist = obj.haversine(grid_pts(:,1), grid_pts(:,2), ...
                                    sat_lat, sat_lon, R_earth);
        
            current_access = dist <= half_swath;
        
            % Detect rises (0 -> 1)'
            rises = current_access & ~data.last_access;
        
            % If this is not first rise, compute revisit
            valid = rises & ~isnan(data.last_rise_time);
            data.revisit_sum(valid) = data.revisit_sum(valid) + ...
                                    (t_current - data.last_rise_time(valid));
            data.revisit_count(valid) = data.revisit_count(valid) + 1;
        
            % Update rise times
            data.last_rise_time(rises) = t_current;
        
            % Update state
            data.last_access = current_access;
        
        
        end
        function get_mean_pass_duration_increment(obj, t_start, t_current, rr, lat_grid, lon_grid, min_elevation, data)
            %%%%%%%Author: Kolja Westphal, ALL RIGHTS RESERVED%%%%%%%%%%%%

            %%%% Input
            % t_start [1x1] start time of the simulation [datetime]
            % t_current [1x1] current time in the simulation [s]
            % rr [3x1] satellite position [km]
            % lat_grid [mxn] latitude grid points [deg]
            % lon_grid [mxn] longitude grid points [deg]
            % half_swath [scalar] half swath width [km]
            % data [PassDuration] Handle class for storing intermediate results
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Constants
            R_earth = 6371;             % km

            % Create a Global Grid (Approx 2-degree spacing) and start NorthWest
            %grid_pts = [deg2rad(lat_grid(:)), deg2rad(lon_grid(:))];
            grid_pts = [deg2rad(lat_grid), deg2rad(lon_grid)];

            % Get current LLA
            t_current_seconds = seconds(t_current);
            time = t_start + t_current_seconds;
            sat_lla = eci2lla(rr'*1000, datevec(time));
            sat_lat = deg2rad(sat_lla(1));
            sat_lon = deg2rad(sat_lla(2));

            % Vectorized Access Calculation
            [dist,c] = obj.haversine(grid_pts(:,1), grid_pts(:,2), ...
                                    sat_lat, sat_lon, R_earth);
        
            % Elevation angle
            e = obj.get_elevation_angle(R_earth, rr, c);
            current_access = min_elevation <= e;
        
            % Detect rises and sets (0 -> 1)'
            rises = current_access & ~data.last_access;
            sets = ~current_access & data.last_access;
        
            % If there was no previous rise, ignore the set
            valid = sets & ~isnan(data.last_rise_time);
            data.pass_sum(valid) = data.pass_sum(valid) + ...
                                    (t_current - data.last_rise_time(valid));
            data.pass_count(valid) = data.pass_count(valid) + 1;
        
            % Update rise times
            data.last_rise_time(rises) = t_current;
        
            % Update state
            data.last_access = current_access;
        
        
        end
    end
end