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
    RequiredEbByNo              = gs_struct.EbByN0;

    %% Frequency Properties
    c              = physconst("LightSpeed");
    lamda          = c / frequency;

    XBAND_UPPER    = 7.75e9;
    XBAND_LOWER    = 7.25e9;
    SBAND_UPPER    = 2.40e9;
    SBAND_LOWER    = 2.20e9;

    if frequency >= SBAND_LOWER && frequency <= SBAND_UPPER
        freqBand = "SBAND";
    elseif frequency >= XBAND_LOWER && frequency <= XBAND_UPPER
        freqBand = "XBAND";
    else
        error("Frequency outside supported bands");
    end

    %% Elevation sweep
    MIN_ELEV = 5;
    MAX_ELEV = 90;
    STEP     = 1;

    angles = MIN_ELEV:STEP:MAX_ELEV;

    %% Transmitted Power
    SAT_EIRP = sat_antenna_gain + sat_transmit_power - sat_cable_loss;

    %% Atmospheric Loss
    L_ion   = 0.1;
    L_atmos = 0.3;

    if freqBand == "XBAND"
        L_rain = 1.1;
    else
        L_rain = 0.3;
    end

    L_env = L_atmos + L_rain + L_ion;

    %% Pointing Loss
    theta       = gs_struct.pointing_loss;
    theta_3dB   = 70 * (lamda / gs_dish_diameter_m);
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
            FSPL    = (lamda/(4*pi*dist))^2;
            FSPL_db = 10*log10(FSPL);

            % Link equation
            Eb_N0_DOWNLINK = SAT_EIRP ...
                + FSPL_db ...
                - L_env ...
                + G_Tsys_Ground ...
                - boltzmann ...
                - dBHz ...
                - L_dp ...
                - gs_cable_loss;

            margins_downlink(i) = Eb_N0_DOWNLINK - RequiredEbByNo;

        end

        % Plot for this orbital height
        plot(angles, margins_downlink, 'LineWidth',2)

        legend_entries(h) = string(orbitalHeight_m/1000) + " km";

    end

    %% Formatting
    % add 3 dB Line
    yline(3, '--', '3 dB limit', 'LineWidth', 2, 'Color', 'r');
    legend(legend_entries, 'Location','best', 'FontSize', 14)
    tit = "Link Budget Margin vs Elevation for Multiple Orbital Heights at " ...
             + "downlink datarate of " + string(datarate / 1e3) + " KBit/s";
    title(tit, "FontSize", 24)
    xlabel("Elevation (deg)","FontSize", 18)
    ylabel("Link Margin (dB)","FontSize", 18)

end