% plotPopulationNoseTrajectories
%

winlen = length(perMouseData(1).traj_window);
ctl_traj = []; ctl_dir = []; run_means_ctl = [];
occr_traj = []; occr_dir = []; run_means_occr = [];
for ii = 3
    tnums = ctl_trials{ii};
    for jj=1:length(tnums)
        ctl_traj = cat(1, ctl_traj, perMouseData(ii).nose_trajectories{tnums(jj)});
        ctl_dir = cat(1, ctl_dir, perMouseData(ii).traj_dir{tnums(jj)});
    end
    tnums = occr_trials{ii};
    for jj=1:length(tnums)
        occr_traj = cat(1, occr_traj, perMouseData(ii).nose_trajectories{tnums(jj)});
        occr_dir = cat(1, occr_dir, perMouseData(ii).traj_dir{tnums(jj)});
    end
end
traj_window = perMouseData(1).traj_window;

%plot some stuff
ltr = ctl_dir > 0;
mean_ltr_ctl = nanmean(ctl_traj(ltr,:));
std_ltr_ctl = nanstd(ctl_traj(ltr,:));
mean_rtl_ctl = nanmean(ctl_traj(~ltr,:));
std_rtl_ctl = nanstd(ctl_traj(~ltr,:));

ltr = occr_dir > 0;
mean_ltr_occr = nanmean(occr_traj(ltr,:));
std_ltr_occr = nanstd(occr_traj(ltr,:));
mean_rtl_occr = nanmean(occr_traj(~ltr,:));
std_rtl_occr = nanstd(occr_traj(~ltr,:));

colors = {'m', [1 .5 1], 'b', [.3 .3 1]};
color_order = [ 1 0 1; [1 .5 1]; 0 0 1; .5 .5 1];
figure;
main_ah = axes('Position', [.1 .1 .6 .8]);
plot(traj_window, mean_ltr_ctl, '-', 'Color', colors{1}, 'LineWidth', 2); hold on;
plot(traj_window, mean_ltr_occr, '--', 'Color', colors{2}, 'LineWidth', 2); 
plot(traj_window, mean_rtl_ctl, '-', 'Color', colors{3}, 'LineWidth', 2);
plot(traj_window, mean_rtl_occr, '--', 'Color', colors{4}, 'LineWidth', 2);
legend_str = {'Control L->R', 'R occluded L->R', 'Control R->L', 'R occluded R->L'};
legend(legend_str, 'Location', 'East');
xlabel('Number of samples @ 40fps'); ylabel('Distance from Trail (px, Rightward is Positive)');
title('Mean Trail Crossing Trajectories');
set(gca, 'YDir', 'reverse');
yl = get(gca, 'ylim');
side_ah = axes('Position', [.75 .1 .15, .8], 'colorOrder', color_order); hold on;
%plot(gca, 1, mean(mean_ltr_ctl), 2, mean(mean_ltr_occr), 3, mean(mean_rtl_ctl), 4, mean(mean_rtl_occr)); 
scatter(gca, 1:4, [mean(mean_ltr_ctl), mean(mean_ltr_occr), mean(mean_rtl_ctl), mean(mean_rtl_occr)], 100, color_order, 'filled'); 
set(gca, 'YDir', 'reverse', 'ylim', yl);
ylabel('Mean Distance from Trail (px)');
