% plot the following segments based on the following information 
[dists, frames] = mt.distanceOnTrail([],1,10);
mt.plotDirection; hold on;
for i=1:size(frames,1)
    range = frames(i,1):frames(i,2);
    np = mt.nosePos(range, :);
    plot(np(:,1), np(:,2), '.m', 'MarkerSize',10);
    %text(mean(np(:,1)), mean(np(:,2)), num2str(dists(i)), 'Color', 'w', 'FontSize', 12);
end