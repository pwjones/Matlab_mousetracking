function ah = popNoseTraj(mouseData, trials,  colorlist, main_ah)
% function ah = popNoseTraj(mouseData)
% 
% Given a data structure from multiple trials (mouseData) this function will plot the 
% averaged nose trajectory during following episodes, with all of the curves on the same
% axis. 

%Setting some defaults for possibly empty arguments
if isempty(main_ah)
    figure;
    main_ah = axes('Position', [.1 .1 .6 .8]);
end
if isempty(trials)
    for ii = 1:length(mouseData)
        trials{ii} = mouseData(ii).file_names;
    end
end
if isempty(colorlist)
    colorlist = {[0 0 0], [0 0 0], [1 0 0],  [1 0 0], [0 0 1], [0 0 1], [0 1 0], [0 1 0]};
end
% For whatever the colors are, the shading is lighter version of each
color_order = [ 0 0 0; 0 0 0; 1 0 0; 1 0 0; 0 0 1; 0 0 1; 0 1 0; 0 1 0];
for ii=1:length(colorlist)
    shading_colorlist{ii} = (colorlist{ii} + [1 1 1])/2;
end

for ii = 1:length(mouseData)
    winlen = length(mouseData(ii).traj_window);
    traj = []; dir = []; run_means = [];

    tnums = trials{ii};
    for jj=1:length(tnums)
        traj = cat(1, traj, mouseData(ii).nose_trajectories{tnums(jj)});
        dir = cat(1, dir, mouseData(ii).traj_dir{tnums(jj)});
    end
    

    %plot some stuff
    ltr = dir > 0;
    mean_ltr = nanmean(traj(ltr,:));    
    std_ltr = nanstd(traj(ltr,:));
    mean_rtl = nanmean(traj(~ltr,:));    
    std_rtl = nanstd(traj(~ltr,:));
    n_ltr = sum(ltr);
    n_rtl = sum(~ltr);

    %figure out the mouse ID
    mouseID = '';
    if ~isempty(mouseData(ii).file_names)
        rem = mouseData(ii).file_names(1); s = '';
        while true
            [s, rem] = strtok(rem, filesep);
            if isempty(rem{1}), break; end
        end
        mouseID = strtok(s,'_');
    end
      
    % Now do the plotting in the main window 
    traj_window = mouseData(ii).traj_window;
    plot_err_poly(main_ah, traj_window, mean_ltr, std_ltr./sqrt(n_ltr), colorlist{ii}, (colorlist{ii}+[1 1 1])/2, 1); hold on;
    plot_err_poly(main_ah, traj_window, mean_rtl, std_rtl./sqrt(n_rtl), colorlist{ii}, (colorlist{ii}+[1 1 1])/2, 1); hold on;
%     plot_err_poly(main_ah, traj_window, mean_ltr_occr, std_ltr_occr./sqrt(n_ltr), colorlist{3}, (colorlist{3}+[1 1 1])/2, 1); hold on;
%     plot_err_poly(main_ah, traj_window, mean_rtl_occr, std_rtl_occr./sqrt(n_rtl), colorlist{4}, (colorlist{4}+[1 1 1])/2, 1); hold on;
%     plot_err_poly(main_ah, traj_window, mean_ltr_occl, std_ltr_occl./sqrt(n_ltr), colorlist{5}, (colorlist{5}+[1 1 1])/2, 1); hold on;
%     plot_err_poly(main_ah, traj_window, mean_rtl_occl, std_rtl_occl./sqrt(n_rtl), colorlist{6}, (colorlist{6}+[1 1 1])/2, 1); hold on;

    %legend_str = {'','Control L->R', '','Control R->L','','R occluded L->R','', 'R occluded R->L','L occluded L->R','', 'L occluded R->L'};
    %legend(legend_str, 'Location', 'East');
    plot(main_ah, traj_window, zeros(size(traj_window)),'--k');
    xlabel('Number of video samples'); ylabel('Distance from Trail (px, Rightward is Positive)');
    title(['Mean Trail Crossing Trajectories - ' mouseID]);
    set(gca, 'YDir', 'reverse', 'TickDir', 'out','FontSize', 14);
    %yl = get(gca, 'ylim');
    %side_ah = axes('Position', [.75 .1 .15, .8], 'colorOrder', color_order, 'TickDir', 'out', 'FontSize', 14); hold on;
    %plot(gca, 1, mean(mean_ltr_ctl), 2, mean(mean_ltr_occr), 3, mean(mean_rtl_ctl), 4, mean(mean_rtl_occr)); 
    %scatter(gca, [1 3 2 4 5 6], [mean(mean_ltr_ctl), mean(mean_rtl_ctl), mean(mean_ltr_occr), mean(mean_rtl_occr), ...
    %                             mean(mean_ltr_occl), mean(mean_rtl_occl)], 100, color_order, 'filled'); 
    %set(gca, 'YDir', 'reverse', 'ylim', yl);
    %ylabel('Mean Distance from Trail (px) - ');

end
