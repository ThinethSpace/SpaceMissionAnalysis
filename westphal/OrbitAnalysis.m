classdef OrbitAnalysis
    methods (Static)
        function distance = haversine(lat1, lon1, lat2, lon2, R)
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
end

    methods
        function mrt = get_mean_revisit_time0(obj, lat_grid, lon_grid, tt, LLA_sat, dt, half_swath)
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
            grid_pts = [deg2rad(lat_grid(:)), deg2rad(lon_grid(:))];
            %grid_pts = [lat_grid, lon_grid];
            num_pts = size(grid_pts, 1);

            sat_lat = deg2rad(LLA_sat(:,1));
            sat_lon = deg2rad(LLA_sat(:,2));

            % Vectorized Access Calculation
            % Pre-allocate a logical matrix: Rows = Grid Points, Cols = Time
            access_matrix = false(num_pts, length(tt));

            for i = 1:length(tt)
                % Haversine distance from satellite to ALL grid points at time i
                dist = obj.haversine(grid_pts(:,1), grid_pts(:,2), sat_lat(i), sat_lon(i), R_earth);
                % Mark grid points within swath
                access_matrix(:, i) = dist <= half_swath;
            end

            mrt_results = zeros(num_pts, 1);

            for p = 1:num_pts
                % Get the visibility timeline for this specific point
                timeline = access_matrix(p, :);
    
                % Find transitions (0 to 1 is a rise, 1 to 0 is a set)
                diffs = diff([0, timeline, 0]);
                rises = find(diffs == 1);
                sets = find(diffs == -1);
    
                if length(rises) > 1
                    % Revisit time is the time from the START of one pass 
                    % to the START of the next pass.
                    revisits = diff(rises) * dt;
                    mrt_results(p) = mean(revisits);
                else
                    mrt_results(p) = NaN; % Not enough passes to calculate mean
                end
            end


            % Reshape for plotting and inverting because grid_points go from west to east
            mrt = reshape(mrt_results, size(lat_grid))';
            %mrt = mrt_results;           

        end
    end
end