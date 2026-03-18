function [] = linkBudget(sat_struct, gs_struct, orbitalHeight_m_arr)

    %% Satellite Properties
    sat_eff_factor              = sat_struct.eff_factor;
    sat_ideal_antenna_gain_db   = sat_struct.ant_gain_ideal;
    sat_antenna_gain            = 10*log10(sat_eff_factor) + sat_ideal_antenna_gain_db;
    sat_transmit_power          = 10*log10(sat_struct.tx_power);
    sat_cable_loss              = sat_struct.loss_cables;
    frequency                   = sat_struct.frequency;
    datarate                    = sat_struct.datarate;

    %% Ground Station Properties
    gs_dish_diameter_m          = gs_struct.dish_diameter;
    gs_altitude_m               = gs_struct.altitude;
    gs_cable_loss               = gs_struct.loss_cables;
    gs_min_elevation            = gs_struct.min_elevation;
    RequiredEbByNo              = gs_struct.EbByN0;


    %% Frequency Properties
    c              = physconst("LightSpeed");
    lambda          = c / frequency;

    XBAND_UPPER    = 7.75e9;
    XBAND_LOWER    = 7.25e9;
    SBAND_UPPER    = 2.40e9;
    SBAND_LOWER    = 2.20e9;

    if frequency >= SBAND_LOWER && frequency <= SBAND_UPPER
        freqBand = "S-Band";
    elseif frequency >= XBAND_LOWER && frequency <= XBAND_UPPER
        freqBand = "X-Band";
    else
        error("Frequency outside supported bands");
    end

    %% Elevation sweep
    MIN_ELEV = gs_min_elevation;
    MAX_ELEV = 90;
    STEP     = 1;

    angles = MIN_ELEV:STEP:MAX_ELEV;

    %% Transmitted Power
    SAT_EIRP = sat_antenna_gain + sat_transmit_power - sat_cable_loss;

    %% Atmospheric Loss
    L_ion   = 0.1;
    L_atmos = 0.3;

    if freqBand == "X-Band"
        L_rain = 1.1;
    else
        L_rain = 0.3;
    end

    L_env = L_atmos + L_rain + L_ion;

    %% Pointing Loss
    theta       = gs_struct.pointing_loss;
    theta_3dB   = 70 * (lambda / gs_dish_diameter_m);
    L_dp        = 12 * (theta/theta_3dB)^2;

    %% Noise
    G_Tsys_Ground = gs_struct.gain_over_temp;
    boltzmann     = 10 * log10(physconst('Boltzmann'));

    %% Data Rate
    dBHz = 10 * log10(datarate);

    %% Plot Setup
    figure
    hold on
    grid on

    %% Loop over orbital heights
    for h = 1:length(orbitalHeight_m_arr)

        orbitalHeight_m = orbitalHeight_m_arr(h);

        margins_downlink = zeros(1,length(angles));

        for i = 1:length(angles)

            % Slant range
            dist = slantRangeCircularOrbit(angles(i), orbitalHeight_m, gs_altitude_m);

            % FSPL
            FSPL    = (lambda/(4*pi*dist))^2;
            FSPL_db =  -10*log10(FSPL); %  Value is negative, adding a minus will make it a positive loss

            % Link equation
            Eb_N0_DOWNLINK = SAT_EIRP ...
                - FSPL_db ...
                - L_env ...
                + G_Tsys_Ground ...
                - boltzmann ...
                - dBHz ...
                - L_dp ...
                - gs_cable_loss;

            margins_downlink(i) = Eb_N0_DOWNLINK - RequiredEbByNo;

        end


        title("Loss Contribution Percentage")

        % Plot for this orbital height
        plot(angles, margins_downlink, 'LineWidth',2)

        legend_entries(h) = string(orbitalHeight_m/1000) + " km";

    end

    % 3 dB requirement line
    yline(3,'--','3 dB requirement',...
        'LineWidth',1.5,...
        'Color',[0.8 0 0],...
        'FontSize',10);

    % Axis labels
    xlabel('Elevation Angle (deg)','FontSize',12)
    ylabel('Link Margin (dB)','FontSize',12)

    % Title
    title(freqBand + " Link Margin vs Elevation",'FontSize',12)

    % Legend
    legend(legend_entries,...
        'Location','northwest',...
        'FontSize',10,...
        'Box','off')

    % Data rate annotation
    text(0.98,0.02,...
        "Data rate: " + string(datarate/1e6) + " Mbit/s",...
        'Units','normalized',...
        'HorizontalAlignment','right',...
        'VerticalAlignment','bottom',...
        'FontSize',10,...
        'BackgroundColor','w',...
        'EdgeColor',[0.7 0.7 0.7])

    % Grid styling
    grid on
    set(gca,...
        'FontSize',11,...
        'LineWidth',1,...
        'GridAlpha',0.3)    %% Formatting

    %% --------------------------------------------------
    %% PROFESSIONAL WATERFALL LINK BUDGET DIAGRAM
    %% (G/T - k grouped to avoid Boltzmann confusion)
    %% --------------------------------------------------

    % Worst case geometry
    angle_loss = MIN_ELEV;
    orbitalHeight_m = orbitalHeight_m_arr(1);

    dist = slantRangeCircularOrbit(angle_loss, orbitalHeight_m, gs_altitude_m);

    % Free Space Path Loss
    FSPL    = (lambda/(4*pi*dist))^2;
    FSPL_db = -10*log10(FSPL); % positive loss

    % Receiver sensitivity term
    GT_minus_k = G_Tsys_Ground - boltzmann;

    steps = [
        SAT_EIRP
        -FSPL_db
        -L_atmos
        -L_rain
        -L_ion
        -L_dp
        -gs_cable_loss
        GT_minus_k
        -dBHz
        -RequiredEbByNo
    ];

    labels = [
        "EIRP"
        "FSPL"
        "Atmosphere"
        "Rain"
        "Ionosphere"
        "Pointing"
        "GS Cable"
        "Receiver (G/T - k)"
        "Data Rate"
        "Required Eb/N0"
    ];

    cumulative = cumsum(steps);

    figure
    hold on
    grid on

    for i = 1:length(steps)

        if steps(i) >= 0
            color = [0.2 0.7 0.2]; % gain
        else
            color = [0.9 0.3 0.3]; % loss
        end

        bar(i, steps(i),'FaceColor',color)

    end

    plot(cumulative,'k-o','LineWidth',2)

    xticks(1:length(labels))
    xticklabels(labels)
    xtickangle(45)

    ylabel("dB")
    title(freqBand + ": Downlink Link Budget Waterfall Diagram")

    margin = cumulative(end);

    plot(length(labels), margin,'kp','MarkerSize',12,'MarkerFaceColor','yellow')

    text(length(labels), margin, sprintf("  Margin = %.2f dB", margin), ...
        'FontSize',12,'VerticalAlignment','bottom')

    ylim([min(cumulative)-10 max(cumulative)+10])

    hold off  
end