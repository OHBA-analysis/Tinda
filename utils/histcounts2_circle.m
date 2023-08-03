function [xCenters, yCenters, counts] = histcounts2_circle(x, y, numBins)

if ~exist('numBins', 'var'), numBins=20; end

range = max(abs([x; y]))*[-1 1];

xEdges = linspace(range(1), range(2), numBins);
yEdges = linspace(range(1), range(2), numBins);

% Calculate the 2D histogram / density map
[counts, xCenters, yCenters] = histcounts2(x, y, xEdges, yEdges);