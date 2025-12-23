% Test: Create porkchop plot

% Create class instances
IT = InterplanetaryTransfers;

% parse data files
[earth_data] = IT.parse_horizon_file('./data/earth_data.txt');
[mars_data] = IT.parse_horizon_file('./data/mars_data.txt');

