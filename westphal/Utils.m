classdef Utils
    methods (Static)

        function plot_orbit_3D(R,Re,plot_title, x_label, y_label, z_label, type)

            
            figure('Name','3D Orbit');

            if type == "dots"
                plot3(R(:,1), R(:,2), R(:,3),'.', 'MarkerSize', 30);
            else
                plot3(R(:,1), R(:,2), R(:,3), 'LineWidth', 1.0);
            end
            hold on; grid on; axis equal;
            % Draw a translucent Earth
            [fX,fY,fZ] = sphere(80);
            surf(Re*fX, Re*fY, Re*fZ, 'FaceAlpha', 0.9, 'EdgeColor', 'none'); colormap gray;

            % Make axes equal and set same range
            axis equal
            all_coords = [R(:,1); R(:,2); R(:,3)];       % combine all data
            range = [min(all_coords) * 1.3, max(all_coords) * 1.3];  % global min/max
            xlim(range); ylim(range); zlim(range);
            constantplane("z",0,FaceAlpha=0.3);

            % Define axis length
            % X-axis
            plot3([0 range(2)], [0 0], [0 0], 'r', 'LineWidth', 1.5);
            % Y-axis
            plot3([0 0], [0 range(2)], [0 0], 'b', 'LineWidth', 1.5);
            % Z-axis
            plot3([0 0], [0 0], [0 range(2)], 'g', 'LineWidth', 1.5);




            %xlabel('x_{ECI} [km]'); ylabel('y_{ECI} [km]'); zlabel('z_{ECI} [km]');
            %title('Unperturbed Keplerian Orbit (ECI/ICRF)');
            xlabel(x_label); ylabel(y_label); zlabel(z_label);
            title(plot_title);
            
            view(35,25);
                 


        end

        function plot_ground_track(lat, lon)

            % Load coastlines
            load coastlines         

            figure
            hold on         

            % Plot land
            plot(coastlon, coastlat, 'w')  % black coastlines           

            % Wrap longitudes to [-180,180] to avoid jumps
            lon = wrapTo180(lon);           

            % Handle dateline jumps by splitting
            jumpIdx = find(abs(diff(lon)) > 180);
            idx = [0; jumpIdx; numel(lon)];         

            % Plot each segment of the ground track
            for k = 1:length(idx)-1
                seg = idx(k)+1 : idx(k+1);
                plot(lon(seg), lat(seg), 'r', 'LineWidth', 1.5)
            end         

            xlabel('Longitude (deg)')
            ylabel('Latitude (deg)')
            xlim([-180 180])
            ylim([-90 90])

            hold off

            

        end
        
        function plot_colormap_earth(lat_range, lon_range, data, title_plot)
            figure;
            % Display the data with NaNs transparent
            h = imagesc(lon_range, lat_range, data);
            set(h, 'AlphaData', ~isnan(data));
                    
            set(gca, 'YDir', 'normal');   % North up
            colormap(jet);
            colorbar;
                    
            title(title_plot);
            xlabel('Longitude'); 
            ylabel('Latitude');
        end

        function plot_heatmap_earth(lat, lon, data, title_plot)

            lat_ranges = -90:2:90;
            lon_ranges = -180:2:180;

            % --- 3. Compute Binned Average ---
            % This bins the data and calculates the sum of Z in each cell
            [N, x_edges, y_edges, x_bin, y_bin] = histcounts2(lon, lat, lon_ranges, lat_ranges);

            res = accumarray([y_bin, x_bin], data, [length(lat_ranges)-1 length(lon_ranges)-1], @(x)mean(x,'omitnan'), NaN);

            figure;
            % Display the data with NaNs transparent
            h = imagesc([-180 180], [-90 90], res);
            set(h, 'AlphaData', ~isnan(res));
                    
            set(gca, 'YDir', 'normal');   % North up
            colormap(jet);
            colorbar;
                    
            title(title_plot);
            xlabel('Longitude'); 
            ylabel('Latitude');

        end

        function plot_ground_station_range(ground_stations, height_apogee, height_perigee, min_elevation, R_earth, additional)

            if ~additional
                figure;
            end
            %psi_max = acos((R_earth/(R_earth+height_apogee))*cos(min_elevation));
            %psi_min = acos((R_earth/(R_earth+height_perigee))*cos(min_elevation));
            %alpha_max = pi/2 - min_elevation - psi_max;
            %alpha_min = pi/2 - min_elevation - psi_min;

            alpha_max = acos((R_earth/(R_earth+height_apogee))*cos(min_elevation)) - min_elevation;
            alpha_min = acos((R_earth/(R_earth+height_perigee))*cos(min_elevation)) - min_elevation;
            hold on;

            for n = 1:size(ground_stations,1)
                for m = 1:2
                    if m == 1
                        psi = alpha_max;
                    else
                        psi = alpha_min;
                    end
                    theta = linspace(0,2*pi,500);

                    lat0 = ground_stations(n,1);   % station latitude
                    lon0 = ground_stations(n,2);  % station longitude

                    lat_circle = asin( sin(lat0)*cos(psi) + ...
                                        cos(lat0)*sin(psi).*cos(theta) );

                    lon_circle = lon0 + atan2( sin(theta).*sin(psi).*cos(lat0), ...
                                                cos(psi)-sin(lat0).*sin(lat_circle) );

                    hold on
                    plot(rad2deg(lon_circle), rad2deg(lat_circle), 'g', 'LineWidth', 2)
                end
            end



        end

        function plot_heatmap_passes(csv)

            csv_file = readmatrix(csv); 
            data = csv_file(2:end,2:end);
            altitudes = csv_file(2:end,1);
            inclinations = csv_file(1,2:end);

            xvalues = string(round(rad2deg(inclinations)));
            yvalues = string(altitudes);
            h = heatmap(xvalues,yvalues,data);


            h.Title = 'Total Pass Duration in hours';
            h.XLabel = 'Inclination (deg)';
            h.YLabel = 'Altitude (km)';
            h.Colormap = nebula;


        end

            
    
    end
end