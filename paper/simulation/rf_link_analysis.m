
%{
    This is a script to analyze the relation between 
    the datarate and orbital height in relation to 
    the link budget
%}




%% XBAND: Satellite Transmitter
sat_xband.eff_factor       = 0.8;               % percentage
sat_xband.ant_gain_ideal   = 14;                % dBi
sat_xband.tx_power         = 2;                 % Watt
sat_xband.loss_cables      = 0.3;               % dB
sat_xband.frequency        = 7.25e9;            % Hz     XBAND
sat_xband.datarate         = 200e6;             % Bps    FINE TUNE THIS ONE


%%  XBAND: Ground Station Properties
% Take the worst performing XBAND ground station from the ESTRACK Network (New Norcia-2)
gs_xband.dish_diameter  = 4.5;        % m
gs_xband.altitude       = 262.94;     % m
gs_xband.gain_over_temp = 28;       % dB / K
gs_xband.pointing_loss  = 0.5;        % deg
gs_xband.loss_cables    = 3;        % dB 
gs_xband.min_elevation  = 20;       % deg
gs_xband.EbByN0         = 5;        % dB QPSK Convolutional 4


%% Analysis

% Altiude and inclination ranges
alt = linspace(180, 400, 5) * 1000; % km

RflinkBudget(sat_xband, gs_xband, alt);