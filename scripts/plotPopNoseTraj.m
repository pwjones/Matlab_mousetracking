% plotPopulationNoseTrajectories
%

winlen = length(perMouseData(1).traj_window);
ctl_traj = []; ctl_dir = []; run_means_ctl = [];
occr_traj = []; occr_dir = []; run_means_occr = [];
occl_traj = []; occl_dir = []; run_means_occl = [];
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
    tnums = occl_trials{ii};
    for jj=1:length(tnums)
        occl_traj = cat(1, occl_traj, perMouseData(ii).nose_trajectories{tnums(jj)});
        occl_dir = cat(1, occl_dir, perMouseData(ii).traj_dir{tnums(jj)});
    end
end
traj_window = perMouseData(1).traj_window;

%plot some stuff
ltr = ctl_dir > 0;
mean_ltr_ctl = nanmean(ctl_traj(ltr,:));
std_ltr_ctl = nanstd(ctl_traj(ltr,:));
mean_rtl_ctl = nanmean(ctl_traj(~ltr,:));
std_rtl_ctl = nanstd(ctl_traj(~ltr,:));
n_ltr = sum(ltr);
n_rtl = sum(~ltr);

ltr = occr_dir > 0;
mean_ltr_occr = nanmean(occr_traj(ltr,:));
std_ltr_occr = nanstd(occr_traj(ltr,:));
mean_rtl_occr = nanmean(occr_traj(~ltr,:));
std_rtl_occr = nanstd(occr_traj(~ltr,:));

ltr = occl_dir > 0;
mean_ltr_occl = nanmean(occl_traj(ltr,:));
std_ltr_occl = nanstd(occl_traj(ltr,:));
mean_rtl_occl = nanmean(occl_traj(~ltr,:));
std_rtl_occl = nanstd(occl_traj(~ltr,:));

plotblue = [.2 .5 1];
%colors = {[1 0 0], [1 0 1], [0 0 1], plotblue};
%color_order = [ 1 0 0; 1 0 1; 0 0 1; .3 .6 1];
colors = {[0 0 0], [0 0 0], [1 0 0], [1 0 0], [0 0 1], [0 0 1]};
color_order = [ 0 0 0; 0 0 0; 1 0 0; 1 0 0; 0 0 1; 0 0 1];

figure;
main_ah = axes('Position', [.1 .1 .6 .8]);
plot_err_poly(main_ah, traj_window, mean_ltr_ctl, std_ltr_ctl./sqrt(n_ltr), colors{1}, (colors{1}+[1 1 1])/2, .5); hold on;
plot_err_poly(main_ah, traj_window, mean_rtl_ctl, std_rtl_ctl./sqrt(n_rtl), colors{2}, (colors{2}+[1 1 1])/2, .5); hold on;
plot_err_poly(main_ah, traj_window, mean_ltr_occr, std_ltr_occr./sqrt(n_ltr), colors{3}, (colors{3}+[1 1 1])/2, .5); hold on;
plot_err_poly(main_ah, traj_window, mean_rtl_occr, std_rtl_occr./sqrt(n_rtl), colors{4}, (colors{4}+[1 1 1])/2, .5); hold on;
plot_err_poly(main_ah, traj_window, mean_ltr_occl, std_ltr_occl./sqrt(n_ltr), colors{5}, (colors{5}+[1 1 1])/2, .5); hold on;
plot_err_poly(main_ah, traj_window, mean_rtl_occl, std_rtl_occl./sqrt(n_rtl), colors{6}, (colors{6}+[1 1 1])/2, .5); hold on;

% plot(traj_window, mean_ltr_ctl, '-', 'Color', colors{1}, 'LineWidth', 2); hold on;
% plot(traj_window, mean_ltr_occr, '--', 'Color', colors{2}, 'LineWidth', 2); 
% plot(traj_window, mean_rtl_ctl, '-', 'Color', colors{3}, 'LineWidth', 2);
% plot(traj_window, mean_rtl_occr, '--', 'Color', colors{4}, 'LineWidth', 2);
legend_str = {'','Control L->R', '','Control R->L','','R occluded L->R','', 'R occluded R->L','L occluded L->R','', 'L occluded R->L'};
legend(legend_str, 'Location', 'East');
plot(main_ah, traj_window, zeros(size(traj_window)),'--k');
xlabel('Number of video samples'); ylabel('Distance from Trail (px, Rightward is Positive)');
title('Mean Trail Crossing Trajectories');
set(gca, 'YDir', 'reverse', 'TickDir', 'out','FontSize', 14);
yl = get(gca, 'ylim');
side_ah = axes('Position', [.75 .1 .15, .8], 'colorOrder', color_order, 'TickDir', 'out', 'FontSize', 14); hold on;
%plot(gca, 1, mean(mean_ltr_ctl), 2, mean(mean_ltr_occr), 3, mean(mean_rtl_ctl), 4, mean(mean_rtl_occr)); 
scatter(gca, [1 3 2 4 5 6], [mean(mean_ltr_ctl), mean(mean_rtl_ctl), mean(mean_ltr_occr), mean(mean_rtl_occr), ...
                             mean(mean_ltr_occl), mean(mean_rtl_occl)], 100, color_order, 'filled'); 
set(gca, 'YDir', 'reverse', 'ylim', yl);
ylabel('Mean Distance from Trail (px)');
