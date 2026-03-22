function plot_stations_dual(gs_rf, gs_opt, rf_names, opt_names)
    % Initialize Figure
    figure('Color', 'w', 'Units', 'inches', 'Position', [1, 1, 9, 6]);
    
    % --- Plot RF Stations (Red Triangles) ---
    h_rf = geoscatter(gs_rf(:,1), gs_rf(:,2), 50, 'red', 'filled', ...
        'Marker', '^', 'MarkerEdgeColor', 'k', 'DisplayName', 'RF Station');
    hold on;
    
    % --- Plot Optical Stations (Blue Circles) ---
    h_opt = geoscatter(gs_opt(:,1), gs_opt(:,2), 50, 'cyan', 'filled', ...
        'Marker', 'o', 'MarkerEdgeColor', 'k', 'DisplayName', 'Optical Station');

    % Set Map Background
    geobasemap grayland % Best for papers (clean and professional)

    % --- Label RF Stations ---
    text(gs_rf(:,1) + 2, gs_rf(:,2), rf_names, 'Color', 'red', ...
        'FontSize', 9, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Interpreter', 'none');

    % --- Label Optical Stations ---
    % Offset these slightly differently to avoid overlap if near an RF station
    text(gs_opt(:,1) - 2, gs_opt(:,2), opt_names, 'Color', [0 0.4 0.7], ...
        'FontSize', 9, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Interpreter', 'none');

    % Add Legend and Title
    legend([h_rf, h_opt], 'Location', 'southoutside', 'Orientation', 'horizontal');
    title('Ground Segment Network: RF and Optical Nodes', 'FontSize', 16);
    
    hold off;
end


%% RF Ground Stations
gs_rf = [37.00 -25.14 275.00; 50.00 5.15 386.68; -31.05 116.19 262.94; ...
         -31.05 116.19 252.26; -35.78 -69.40 1550.00; 5.25 -52.80 -14.67; ...
         67.86 20.97 400.68; 67.86 20.96 402.17; 40.45 -4.37 794.09];

gs_rf_names = {'Santa Maria'; 'Redu'; 'New Norcia-2'; 'New Norcia-1'; ...
         'Malargue'; 'Kourou'; 'Kiruna-2'; 'Kiruna-1'; 'Cebreros'};

%% Optical Ground Stations
gs_opt = [
  37.84  22.62  440
  35.33  25.13  400
  28.30  -16.5 2400
  37.09  -2.36  489
  52.85  10.13   72
  48.08  11.64  549
  53.33  13.07   66
  48.08  11.28  615
];

opt_names = {'Namea', 'Heraklion', 'Teide', 'FOGATA', 'LaBoT', ...
             'OGS NBB', 'OGS NSG', 'OGSOP-NG'};

%% Run the plot
plot_stations_dual(gs_rf, gs_opt, gs_rf_names, opt_names);