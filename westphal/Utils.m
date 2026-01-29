classdef Utils
    methods (Static)

        function plot_orbit_3D(R,Re,title, x_label, y_label, z_label)

            figure('Name','3D Orbit');

            plot3(R(:,1), R(:,2), R(:,3), 'LineWidth', 1.0); hold on; grid on; axis equal;
            % Draw a translucent Earth
            [fX,fY,fZ] = sphere(40);
            surf(Re*fX, Re*fY, Re*fZ, 'FaceColor','w', 'EdgeColor','none');

            % Coastline
            load coastlines
            lat = deg2rad(coastlat);
            lon = deg2rad(coastlon);

            xc = Re*cos(lat).*cos(lon);
            yc = Re*cos(lat).*sin(lon);
            zc = Re*sin(lat);

            % Natural green coast
            plot3(xc, yc, zc, 'k', 'LineWidth', 1.5)

            %xlabel('x_{ECI} [km]'); ylabel('y_{ECI} [km]'); zlabel('z_{ECI} [km]');
            %title('Unperturbed Keplerian Orbit (ECI/ICRF)');
            xlabel(x_label); ylabel(y_label); zlabel(z_label);
            title(title);
            
            view(35,25);

        end

        function plot_ground_track(lat, lon)

            figure
            geobasemap('darkwater')
            hold on
                    
            % Wrap longitudes
            lon = wrapTo180(lon);
                    
            % Find big jumps (>180°) in longitude
            jumpIdx = find(abs(diff(lon)) > 180);
                    
            % Split indices into segments
            idx = [0; jumpIdx; numel(lon)];
                    
            % Plot each segment separately
            for k = 1:length(idx)-1
                seg = idx(k)+1 : idx(k+1);
                geoplot(lat(seg), lon(seg), 'r', 'LineWidth', 1.5)
            end

            

        end
    end
end