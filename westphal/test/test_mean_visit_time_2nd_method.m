%% 1. Create a meshgrid of latitude and longitude
nLat = 20; nLon = 20;          % number of layers
latVec = linspace(-90, 90, nLat);   % latitude range
lonVec = linspace(-180, 180, nLon); % longitude range

[LonGrid, LatGrid] = meshgrid(lonVec, latVec);

%% 2. Define some example longitude intervals
% Each row: [startLon, endLon]
lonIntervals = [
    -20, 20;
    -60, 60;
    -100, 100
];
latIntervals = [-50, 50; -30, 30; 10, 10];

%% 3. Initialize a matrix to check coverage
Coverage = zeros(size(LonGrid));

%% 4. Check which points are inside any interval
for k = 1:size(lonIntervals,1)
    startLon = lonIntervals(k,1);
    endLon = lonIntervals(k,2);
    startLat = latIntervals(k,1);
    endLat= latIntervals(k,2);

    
    
    % Account for wrap-around at 180/-180 if needed
    inLatInterval = (latVec >= startLat) & (latVec <= endLat);
    inLonInterval = (lonVec >= startLon) & (lonVec <= endLon);

    
    % Mark points covered
    % Expand 1D vectors to 2D grid using ndgrid
    [LatMask, LonMask] = ndgrid(inLatInterval, inLonInterval);

    % Combine and mark Coverage
    Coverage = Coverage | (LatMask & LonMask);  % both lat and lonCoverage = Coverage | inInterval;
end

%% 5. Display results
disp('Coverage matrix (1 = covered, 0 = not covered):');
disp(Coverage);

%% 6. Optional: visualize
figure;
imagesc(lonVec, latVec, Coverage);
xlabel('Longitude'); ylabel('Latitude');
title('Coverage Map');
colorbar;
