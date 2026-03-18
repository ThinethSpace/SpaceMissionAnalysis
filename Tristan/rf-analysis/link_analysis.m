
%{
    This is a script to analyze the relation between 
    the datarate and orbital height in relation to 
    the link budget
%}

%% SBAND: Satellite Transmitter
sat_sband.eff_factor       = 0.8;       % percentage
sat_sband.ant_gain_ideal   = 7;         % dBi
sat_sband.tx_power         = 1;         % Watt
sat_sband.loss_cables      = 0.3;       % dB
sat_sband.frequency        = 2.2e9;     % Hz SBAND
sat_sband.datarate         = 2e6;       % Bps    FINE TUNE THIS ONE


%%  SBAND: Ground Station Properties
% Take the worst performing SBAND ground station from the ESTRACK Network (Santa Maria)
gs_sband.dish_diameter      = 5;        % m
gs_sband.altitude           = 275;      % m
gs_sband.gain_over_temp     = 16;       % dB / K
gs_sband.pointing_loss      = 1;        % deg
gs_sband.loss_cables        = 3;        % dB 
gs_sband.min_elevation      = 10;       % Minimum allowed elevation for downlink
gs_sband.EbByN0             = 5;        % dB QPSK Convolutional 4


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

alt = 200e3:50e3:400e3;

linkBudget(sat_sband, gs_sband, alt);
linkBudget(sat_xband, gs_xband, alt)