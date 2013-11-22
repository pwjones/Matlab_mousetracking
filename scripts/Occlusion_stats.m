% The statistics on if things change during occlusion

%matrix structure - columns are rew/distractor, rows are control/occluded
propFollowed = cell(2,2);
propFollowed{1,1} = catStructArray(perMouseData, 'rew_propFollowed', ctl_trials);
propFollowed{1,2} = catStructArray(perMouseData, 'dist_propFollowed', ctl_trials);
propFollowed{2,1} = catStructArray(perMouseData, 'rew_propFollowed', occ_trials);
propFollowed{2,2} = catStructArray(perMouseData, 'dist_propFollowed', occ_trials);

[p, h] = ranksum(propFollowed{1,1}, propFollowed{2,1});
if h sig = ''; else sig = 'NOT'; end
sig_str = sprintf('Rewarded Unoccluded/Occ Following Proportions Wilcoxon: is %s significantly different: h=%i p=%f', sig, hr, p);
disp(sig_str);

%% Trail Area Explored
%matrix structure - columns are rew/distractor, rows are control/occluded
distFollowed = cell(2,2);
distFollowed{1,1} = catStructArray(perMouseData, 'rew_trail_area', ctl_trials) .* propFollowed{1,1};
distFollowed{1,2} = catStructArray(perMouseData, 'dist_trail_area', ctl_trials) .* propFollowed{1,2};
distFollowed{2,1} = catStructArray(perMouseData, 'rew_trail_area', occ_trials) .* propFollowed{2,1};
distFollowed{2,2} = catStructArray(perMouseData, 'dist_trail_area', occ_trials) .* propFollowed{2,2};

[p, h] = ranksum(distFollowed{1,1}, distFollowed{2,1});
if h sig = ''; else sig = 'NOT'; end
sig_str = sprintf('Rewarded Unoccluded/Occ Following Distances Wilcoxon: is %s significantly different: h=%i p=%f', sig, hr, p);
disp(sig_str);


%% Trail Area Exploration Rate
%matrix structure - columns are rew/distractor, rows are control/occluded
rateFollowed = cell(2,2);
rateFollowed{1,1} = distFollowed{1,1} ./ (catStructArray(perMouseData, 'total_frames', ctl_trials) ./ catStructArray(perMouseData, 'frame_rate', ctl_trials));
rateFollowed{1,2} = distFollowed{1,2} ./ (catStructArray(perMouseData, 'total_frames', ctl_trials) ./ catStructArray(perMouseData, 'frame_rate', ctl_trials));
rateFollowed{2,1} = distFollowed{2,1} ./ (catStructArray(perMouseData, 'total_frames', occ_trials) ./ catStructArray(perMouseData, 'frame_rate', occ_trials));
rateFollowed{2,2} = distFollowed{2,2} ./ (catStructArray(perMouseData, 'total_frames', occ_trials) ./ catStructArray(perMouseData, 'frame_rate', occ_trials));

[p, h] = ranksum(rateFollowed{1,1}, rateFollowed{2,1});
if h sig = ''; else sig = 'NOT'; end
sig_str = sprintf('Rewarded Unoccluded/Occ Following Rates Wilcoxon: is %s significantly different: h=%i p=%f', sig, hr, p);
disp(sig_str);

[p, h] = ranksum(rateFollowed{1,2}, rateFollowed{2,2});
if h sig = ''; else sig = 'NOT'; end
sig_str = sprintf('Distractor Unoccluded/Occ Following Rates Wilcoxon: is %s significantly different: h=%i p=%f', sig, hr, p);
disp(sig_str);


%% Plot the values that we have computed

g = {'Rewarded Unoccluded', 'Rewarded Occluded', 'Distractor Unoccluded', 'Distractor Occluded'};
% Need to find a matrix size to do a box plot
maxTrials = 0;
for ii = 1:size(rateFollowed, 1)
    if length(rateFollowed{ii,1}) > maxTrials;    maxTrials = length(rateFollowed{ii,1});    end
end
% Reformatting into matrix
rateFollowedM = NaN*zeros(maxTrials, 4);
k=1;
for ii = 1:2
    for jj = 1:2
        temp = rateFollowed{ii,jj};
        rateFollowedM(1:length(temp), k) = temp;
        temp = distFollowed{ii,jj};
        distFollowedM(1:length(temp), k) = temp;
        k = k+1;
    end
end

figure;
mean_rates = reshape(nanmean(rateFollowedM)', 2,2);
std_rates = reshape(nanstd(rateFollowedM)' ./ sqrt(size(rateFollowedM,1))', 2,2);
bar(mean_rates);
set(gca, 'TickDir', 'out', 'FontSize', 14);
ylabel('Trail Following Rate, mm^2/sec', 'FontSize', 16);
addErrorBars(gca, [.85 1.15 1.85 2.15], mean_rates', std_rates', 'k', .1);

figure;
mean_dists = reshape(nanmean(distFollowedM)', 2,2);
std_dists = reshape(nanstd(distFollowedM)' ./ sqrt(size(distFollowedM,1))', 2,2);
bar(mean_dists);
set(gca, 'TickDir', 'out', 'FontSize', 14);
ylabel('Trail Following Area, mm^2/sec', 'FontSize', 16);
addErrorBars(gca, [.8 1.2 1.8 2.2], mean_dists', std_dists', 'k', .1);

