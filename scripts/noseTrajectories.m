function [nose_traj, dirs, wind, crossingInds] = noseTrajectories(mt, varargin)
% function [nose_traj, wind, dir] = noseTrajectories(mt)
%
% Want to try to look at the nose trajectories for when the animal is tracking the trail.
% Basically, want to look at the lateral distance from the trail after when the animal sweeps 
% left->right or right->left.
if length(varargin) >= 1 && ~isempty(varargin{1})
    wind = varargin{1};
else
    %wind = -10:20; %window for looking at the positions
    wind = -10:30; %window for looking at the positions
end
if length(varargin) >= 2
    pb = varargin{2};
else
    pb = 0;
end
[crossings, dirs] = mt.findTrailCrossings([],1); %finds the crossings and their directions
crossingInds = crossings;
% let's plot them on the trail for sanity check
if pb
    mt.plotPosition([]);
    np = mt.nosePos(crossings, :);
    plot(np(:,1), np(:,2), 'rx');
end

nose_traj = NaN*zeros(length(crossings), length(wind));
for ii = 1:length(crossings)
    cross_wind = wind+crossings(ii);
    valid = cross_wind > 0 & cross_wind <= mt.nFrames;
    nose_traj(ii,valid) = mt.orthogonalDistFromTrail(cross_wind(valid),1);
end

ltr = dirs > 0;
rtl = dirs < 1;
mean_ltr = nanmean(nose_traj(ltr, :));
mean_rtl = nanmean(nose_traj(rtl, :));

if pb
    figure; hold on;
    plot(wind', nose_traj(ltr,:)', 'm');
    plot(wind', nose_traj(rtl,:)', 'b');
    plot(wind', mean_ltr, 'm', 'LineWidth', 2);
    plot(wind', mean_rtl, 'b', 'LineWidth', 2);
end
