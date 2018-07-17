function [ah, propTurnsLeft]  = makeTurningHistogram(allDirs, allTurnPos, binEdgeVect, varargin)
% function makeTurningHistogram(allDirs, allTurnPos, binEdgeVect)
% allDirs - direction is coded as follows, -1: rightward, 1: leftward
% allTurnPos - position is relative to trail, rightward is positive,
% leftward negative
% binEdgeVect - a vector of position bin edges
% varargin{1} - axis handle for plotting
% varargin{2} - Linespec string for plotting options

%edgeVect = -20:2:20;
[N,bin] = histc(allTurnPos, binEdgeVect);
nbins = length(N);
propTurnsLeft = zeros(nbins, 1);
for ii = 1:nbins
    temp = allDirs(bin == ii); % turns at each distance
    % num leftward turns / num all turns
    propTurnsLeft(ii) = sum(temp == 1) ./ length(temp); 
end

% now just check args and plot
if nargin<4
    figure;
    ah = axes;
else
    ah = varargin{1};
end
if nargin<5
    linespec = 'k-';
else
    linespec = varargin{2};
end

lh = stairs(ah, binEdgeVect, propTurnsLeft, linespec);
set(lh, 'lineWidth', 2);
xlabel('Distance from Trail (mm)');
ylabel('Prop Turns Leftwards');
xlim([min(binEdgeVect) max(binEdgeVect)]);