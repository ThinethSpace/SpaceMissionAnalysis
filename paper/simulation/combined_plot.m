%% Combined RF (X-Band) + Optical Link Margin Plot
% Orbit heights: 200–450 km, elevation: 10–90°, inline labels

clear; clc;

%% ── SHARED GEOMETRY ────────────────────────────────────────────────────
Re          = 6371e3;           % Earth radius (m)
elev_deg    = 10:1:90;
elev_rad    = deg2rad(elev_deg);
heights_km  = 200:50:450;       % [200 250 300 350 400 450] km
heights_m   = heights_km * 1e3;

% Slant range (m) — shared formula
slant_fn = @(elev_r, h_m) ...
    sqrt((Re*sin(elev_r)).^2 + 2*Re*h_m + h_m^2) - Re*sin(elev_r);

%% ── COLOURS ─────────────────────────────────────────────────────────────
% All RF lines: one blue shade
% All Optical lines: one orange shade
% Both are solid; text labels carry the altitude info
rf_color  = [0.15 0.45 0.80];   % steel blue
opt_color = [0.90 0.45 0.05];   % burnt orange

label_elev_deg = 50;            % elevation angle at which to place inline label
[~, label_idx] = min(abs(elev_deg - label_elev_deg));

%% ══════════════════════════════════════════════════════════════════════
%% 1.  RF  –  X-Band Link Margin
%% ══════════════════════════════════════════════════════════════════════

% ── Satellite (X-Band) ────────────────────────────────────────────────
frequency       = 8.0e9;        % Hz  (X-Band centre)
c               = physconst('LightSpeed');
lambda_rf       = c / frequency;

eff_factor      = 0.80;
ant_gain_ideal  = 7;            % dBi
ant_gain_dB     = 10*log10(eff_factor) + ant_gain_ideal;
tx_power_dBW    = 10*log10(1);  % 1 W → 0 dBW
cable_loss_sat  = 1;            % dB

SAT_EIRP        = ant_gain_dB + tx_power_dBW - cable_loss_sat;

% ── Ground Station ────────────────────────────────────────────────────
gs_altitude_m   = 0;
dish_diam_m     = 2.4;          % m
G_T             = 20;           % dB/K  (G/T)
cable_loss_gs   = 1;            % dB
pointing_err    = 0.2;          % deg
theta_3dB       = 70 * (lambda_rf / dish_diam_m);
L_pointing      = 12 * (pointing_err / theta_3dB)^2;  % dB

RequiredEbN0    = 9.6;          % dB  (e.g. QPSK, BER 1e-6)
datarate        = 100e6;        % 100 Mbit/s
dBHz            = 10*log10(datarate);
boltzmann_dB    = 10*log10(physconst('Boltzmann'));

% ── Atmospheric losses (X-Band) ───────────────────────────────────────
L_atmos = 0.3;   % dB
L_rain  = 1.1;   % dB  (X-Band)
L_ion   = 0.1;   % dB
L_env   = L_atmos + L_rain + L_ion;

%% ══════════════════════════════════════════════════════════════════════
%% 2.  OPTICAL  –  1550 nm Link Margin
%% ══════════════════════════════════════════════════════════════════════

lambda_opt   = 1550e-9;         % m

P_tx_W       = 2.5;             % W
P_tx_dBm     = 10*log10(P_tx_W * 1e3);   % dBm

D_tx         = 0.08;            % m  (8 cm satellite aperture)
D_rx         = 0.70;            % m  (70 cm ground aperture)
eta          = 0.75;            % aperture efficiency (FIX: was 1.0)

G_tx_dB = 10*log10(eta) + 20*log10((pi*D_tx) / lambda_opt);
G_rx_dB = 10*log10(eta) + 20*log10((pi*D_rx) / lambda_opt);

L_pointing_opt  = 4;            % dB  pointing loss
L_turbulence    = 4;            % dB  atmospheric scintillation
L_hardware      = 5;            % dB  optical train + detector losses
L_lumped        = L_pointing_opt + L_turbulence + L_hardware; % = 13 dB
                                % (FIX: was 20 dB — too pessimistic as single value)

L_zenith     = 0.46;           % dB  zenith atmospheric loss

% FIX: sensitivity — use -22 dBm for 10 Gbps IM-DD (conservative)
% -35 dBm is coherent and overly optimistic for this scenario
Sensitivity  = -22;            % dBm

%% ══════════════════════════════════════════════════════════════════════
%% 3.  FIGURE
%% ══════════════════════════════════════════════════════════════════════

figure('Color','w','Position',[100 100 1000 620]);
hold on; grid on; box on;

%% ── RF loop ──────────────────────────────────────────────────────────
for k = 1:length(heights_km)
    h_m  = heights_m(k);
    h_km = heights_km(k);

    margin_rf = zeros(1, length(elev_deg));
    for i = 1:length(elev_deg)
        dist   = slant_fn(elev_rad(i), h_m);
        FSPL   = (lambda_rf / (4*pi*dist))^2;
        FSPL_dB = -10*log10(FSPL);

        EbN0 = SAT_EIRP ...
             - FSPL_dB ...
             - L_env ...
             + G_T ...
             - boltzmann_dB ...
             - dBHz ...
             - L_pointing ...
             - cable_loss_gs;

        margin_rf(i) = EbN0 - RequiredEbN0;
    end

    plot(elev_deg, margin_rf, ...
         'Color', rf_color, 'LineWidth', 1.8);

    % Inline label
    y_label = margin_rf(label_idx);
    text(elev_deg(label_idx), y_label, ...
         sprintf(' %d km', h_km), ...
         'Color', rf_color, ...
         'FontSize', 8, ...
         'FontWeight', 'bold', ...
         'VerticalAlignment', 'middle', ...
         'Clipping', 'on');
end

%% ── Optical loop ─────────────────────────────────────────────────────
for k = 1:length(heights_km)
    h_m  = heights_m(k);
    h_km = heights_km(k);

    slant_m  = slant_fn(elev_rad, h_m);        % vector over elevation

    % Free-Space Path Loss (dB)
    L_path_dB = 20*log10(lambda_opt ./ (4*pi .* slant_m));   % negative value = loss

    % Atmospheric loss (zenith scaled by airmass)
    L_atm_dB  = L_zenith ./ sin(elev_rad);

    % Received power (dBm)
    P_rx_dBm = P_tx_dBm + G_tx_dB + L_path_dB ...
               - L_atm_dB + G_rx_dB - L_lumped;

    Margin = P_rx_dBm - Sensitivity;

    plot(elev_deg, Margin, ...
         'Color', opt_color, 'LineWidth', 1.8);

    % Inline label — offset vertically so it doesn't clash with RF label
    y_label = Margin(label_idx);
    text(elev_deg(label_idx), y_label, ...
         sprintf(' %d km', h_km), ...
         'Color', opt_color, ...
         'FontSize', 8, ...
         'FontWeight', 'bold', ...
         'VerticalAlignment', 'middle', ...
         'Clipping', 'on');
end

%% ── 3 dB requirement line ────────────────────────────────────────────
yline(3, '--', 'LineWidth', 1.5, 'Color', [0.6 0 0], ...
      'HandleVisibility', 'off');
text(91, 3, '3 dB req.', 'Color', [0.6 0 0], ...
     'FontSize', 9, 'VerticalAlignment', 'middle', 'Clipping', 'off');

%% ── Dummy lines for legend ───────────────────────────────────────────
h_rf  = plot(nan, nan, '-',  'Color', rf_color,  'LineWidth', 2.5);
h_opt = plot(nan, nan, '-',  'Color', opt_color, 'LineWidth', 2.5);
legend([h_rf, h_opt], ...
       {'X-Band RF (100 Mbit/s, 1W, QPSK)', ...
        'Optical 1550 nm (10 Gbps, 2.5W)'}, ...
       'Location', 'northwest', 'FontSize', 11, 'Box', 'off');

%% ── Axes formatting ──────────────────────────────────────────────────
xlabel('Elevation Angle (°)',  'FontSize', 13);
ylabel('Link Margin (dB)',     'FontSize', 13);
title('Combined RF & Optical Link Margin — VLEO to Ground Station', ...
      'FontSize', 14);
subtitle(sprintf(['X-Band: P_{tx}=1W, f=8GHz  |  ' ...
                  'Optical: P_{tx}=2.5W, \\lambda=1550nm, ' ...
                  'D_{tx}=8cm, D_{rx}=70cm']), ...
         'FontSize', 9, 'Color', [0.4 0.4 0.4]);

xlim([10 90]);
ax = gca;
ax.FontSize  = 11;
ax.LineWidth = 1;
ax.GridAlpha = 0.3;

% Annotation box explaining colours
annotation('textbox', [0.68 0.03 0.22 0.12], ...
    'String', {'\color[rgb]{0.15,0.45,0.80}■ X-Band (200→450 km)', ...
               '\color[rgb]{0.90,0.45,0.05}■ Optical (200→450 km)', ...
               'Labels: altitude per line'}, ...
    'FontSize', 9, 'EdgeColor', [0.8 0.8 0.8], ...
    'BackgroundColor', 'w', 'FitBoxToText', 'on', ...
    'Interpreter', 'tex');