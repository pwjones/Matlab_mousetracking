% Compile the turn averages of the trials in exp

wind = -20:20;
ci = 21;
followingThresh = 20;
dirs = []; dists = [];
for ii = 1:length(exp.vids)
    [tf, dr, di] = exp.vids(ii).findFollowingTurns([], 1, followingThresh, wind);
    dirs = cat(1, dirs, dr);
    dists = cat(2, dists, di);
end

pos = dirs == 1;
peaks = dists(ci, :);
pos_peaks = dists(ci, pos);
pos_peak_mean = mean(squeeze(pos_peaks));
neg = dirs == -1;
neg_peaks = dists(ci, neg);
neg_peak_mean = mean(squeeze(neg_peaks));
recentered_dists = zeros(size(dists));
scat = rand(size(dists, 2), 1)*.1*pos_peak_mean*2 - pos_peak_mean;
for ii = 1:size(dists, 2)
    if (dirs(ii) == 1)
        recentered_dists(:,ii) = (dists(:,ii) - peaks(ii)) + pos_peak_mean;
    else
        recentered_dists(:,ii) = (dists(:,ii) - peaks(ii)) + neg_peak_mean;
    end
end
pi = find(pos);
mean_pos = nanmean(recentered_dists(:,pos),2);
mean_neg = nanmean(recentered_dists(:,neg),2);
xt = wind ./ exp.vids(1).frameRate * 1000;
figure; 
plot(xt, zeros(size(xt)),'--k'); hold on;
plot(xt, recentered_dists(:,1:40:end), 'Color', [.5 .5 .5]); 
plot(xt, mean_pos, 'LineWidth', 2); hold on; plot(xt, mean_neg,'g', 'LineWidth', 2);
xlim([-200 200]);
ylim([2*neg_peak_mean  2*pos_peak_mean]);
set(gca, 'Ydir', 'reverse');
xlabel('Time Relative to Turn (msec)');
ylabel('Distance From Trail (px)');