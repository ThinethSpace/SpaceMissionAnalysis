classdef Utils
    methods (Static)

        function plot_orbit_3D(R,Re,plot_title, x_label, y_label, z_label)

            
            figure('Name','3D Orbit');

            plot3(R(:,1), R(:,2), R(:,3), 'LineWidth', 1.0); 
            hold on; grid on; axis equal;
            % Draw a translucent Earth
            %[fX,fY,fZ] = sphere(80);
            %surf(Re*fX, Re*fY, Re*fZ, 'FaceAlpha', 0.2, 'EdgeColor', 'none'); colormap gray;

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

            

        end
    end
end