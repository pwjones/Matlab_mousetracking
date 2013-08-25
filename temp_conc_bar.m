barMat = [nanmean(rew_prop50), nanmean(rew_prop); nanmean(dist_prop50), nanmean(dist_prop); nanmean(fake_prop50), nanmean(fake_prop)];
bar(barMat);

figure;
all_rew_dists = []; all_distract_dists = [];
for ii = 1:length(rew_dists)
    all_rew_dists = cat(2, all_rew_dists, rew_dists{ii});
    all_distract_dists = cat(2, all_distract_dists, distract_dists{ii});
end
all_rew_dists50 = []; all_distract_dists50 = [];
for ii = 1:length(rew_dists50)
    all_rew_dists50 = cat(2, all_rew_dists50, rew_dists50{ii});
    all_distract_dists50 = cat(2, all_distract_dists50, distract_dists50{ii});
end
barMat = [nanmean(all_rew_dists50), nanmean(all_rew_dists); nanmean(all_distract_dists50), nanmean(all_distract_dists)];
bar(barMat);
