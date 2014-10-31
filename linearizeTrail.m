function linearizeTrail(tracker)

dist_thresh = 50;

ff = tracker.getFollowingSegments([],1, dist_thresh);
nsegs = size(ff,1);
figure; hold on;
dist_ax = subplot(2,1,1); hold on;
ylim([-dist_thresh*1.5 dist_thresh*1.5]);
xlim([1 tracker.nFrames]);
vel_ax = subplot(2,1,2); hold on;
xlim([1 tracker.nFrames]);
plot(dist_ax, 1:tracker.nFrames, zeros(1, tracker.nFrames), 'k', 'LineWidth', 2);
for i=1:nsegs
    frames = ff(i,1):ff(i,2);
    dists = tracker.orthogonalDistFromTrail(frames,1);
    crossings = tracker.findTrailCrossings(frames,1);
    plot(dist_ax, frames, dists, '-','LineWidth', 1, 'Color', 'b');
    plot(dist_ax, frames(crossings), dists(crossings), 'rx', 'MarkerSize', 10);
    vels = tracker.noseVel(frames);
    plot(vel_ax, frames, vels, '-k', 'LineWidth', 2);
end
ylim([-dist_thresh*1.5 dist_thresh*1.5]);

