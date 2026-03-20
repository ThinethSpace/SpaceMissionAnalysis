classdef OrbitAnalysis
    methods (Static)
        function [distance,c] = haversine(lat1, lon1, lat2, lon2, R)
            %%%%%%%Author: Kolja Westphal, ALL RIGHTS RESERVED%%%%%%%%%%%%
            
            %%%% Input
            % lat1, lon1 [nx1] point 1 latitude and longitude [rad]
            % lat2, lon2 [nx1] point 2 latitude and longitude [rad]
            % R [1x1] radius of the sphere 
            
            %%%% Output 
            % distance [1x1] distance between the two points 
            % c [1x1] central angle between the two points [rad]
            
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
            e = atan2(r_sat .* cos(central_angle) - R_earth, r_sat .* sin(central_angle));
        end
end

    methods
        function mrt = get_mean_revisit_time(obj, lat_grid, lon_grid, tt, LLA_sat, half_swath, R_earth)
            %%%%%%%Author: Kolja Westphal, ALL RIGHTS RESERVED%%%%%%%%%%%%
            
            %%%% Input
            % lat_grid [mxn] latitude grid points [deg]
            % lon_grid [mxn] longitude grid points [deg]
            % tt [nx1] time vector [sec]
            % LLA_sat [nx3] satellite lat, lon, alt [deg, deg, km]
            % half_swath [1x1] half swath width [km]

            %%%% Output
            % mrt [1x1] mean revisit time [sec]
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Conert grid to radians
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
        function get_mean_revisit_time_increment(obj, t_start, t_current, rr, lat_grid, lon_grid,half_swath,R_earth, data) 
            %%%%%%%Author: Kolja Westphal, ALL RIGHTS RESERVED%%%%%%%%%%%%
            
            %%%% Input
            % t_start [1x1] start time of the simulation [datetime]
            % t_current [1x1] current time in the simulation [s]
            % rr [3x1] satellite position [km]
            % lat_grid [mxn] latitude grid points [deg]
            % lon_grid [mxn] longitude grid points [deg]
            % half_swath [1x1] half swath width [km]
            % R_earth [1x1] Earth radius [km]
            % data [SharedData] struct for storing intermediate results
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Convert grid to radians
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
        function get_mean_pass_duration_increment(obj, t_start, t_current, rr, lat_grid, lon_grid, min_elevation, R_earth, optical_availability, data)
            %%%%%%%Author: Kolja Westphal, ALL RIGHTS RESERVED%%%%%%%%%%%%

            %%%% Input
            % t_start [1x1] start time of the simulation [datetime]
            % t_current [1x1] current time in the simulation [s]
            % rr [3x1] satellite position [km]
            % lat_grid [mxn] latitude grid points [deg]
            % lon_grid [mxn] longitude grid points [deg]
            % R_earth [1x1] Earth radius [km]
            % data [PassDuration] Handle class for storing intermediate results
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Convert grid to radians
            grid_pts = [deg2rad(lat_grid), deg2rad(lon_grid)];

            % Get current LLA
            t_current_seconds = seconds(t_current);
            time = t_start + t_current_seconds;
            sat_lla = eci2lla(rr'*1000, datevec(time));
            sat_lat = deg2rad(sat_lla(1));
            sat_lon = deg2rad(sat_lla(2));



            % Vectorized Access Calculation
            [~,c] = obj.haversine(grid_pts(:,1), grid_pts(:,2), ...
                                    sat_lat, sat_lon, R_earth);
        
            % Elevation angle
            e = obj.get_elevation_angle(R_earth, rr, c);
            current_access = min_elevation <= e;
        
            % Detect aos and los
            aos = current_access & ~data.last_access;
            los = ~current_access & data.last_access;

            % Check if there has been aos before
            signaling = current_access & data.last_access;

            % Update FOM for currently accessed ground stations
            data.fom(signaling,1) = data.fom(signaling,1) + e(signaling);
            data.fom(signaling,2) = data.fom(signaling,2) + e(signaling) .* (t_current - data.last_aos_time(signaling));
            % add here the fom for weather
        
            % If there was no previous rise, ignore the set
            valid_los = los & ~isnan(data.last_aos_time);
            
            % Evaluate which LOS to accept based on FOM
            accepted_los = false(size(valid_los));
            current_fom = zeros(size(valid_los));
            combined = current_access | valid_los;
            if any(valid_los)
                % Determine the total fom based from the components
                current_fom(combined) = data.fom(combined,1) -  (1 / t_current - data.last_aos_time(combined) .* data.fom(combined,2));

                % Add cloud probability here if available
                % If clouds are present, set the FOM to a negative value to prevent selection of this LOS
                current_fom(combined) = current_fom(combined) .*(2*(rand(size(current_fom(combined))) < optical_availability(combined)) - 1);

                % Find the GS that currently has the highest FOM
                [~, idx] = max(current_fom);


                % If the highest FOM is from a valid LOS, accept it
                % Here, if only one LOS is valid, it can be still negative (due to cloads), so maximum FOM would then be some zero number
                if valid_los(idx) == 1

                    % Accept this LOS and reset FOM for all ground stations and set the last AOS for current accessed GS to current time, simulating a new AOS with that GS
                    accepted_los(idx) = 1; 
                    data.fom(:,:) = 0;
                    data.last_aos_time(current_access) = t_current;
                end
            end

            % Update pass duration for accepted LOS
            data.pass_sum(accepted_los) = data.pass_sum(accepted_los) + ...
                                (t_current - data.last_aos_time(accepted_los));
            data.pass_count(accepted_los) = data.pass_count(accepted_los) + 1;

            % Update AOS times for newly accessed ground stations
            data.last_aos_time(aos) = t_current;
        
            % Update state
            data.last_access = current_access;
        
        end
    end
end