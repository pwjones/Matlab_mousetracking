% popTurningAnal.m

% population analysis of the turning behavior
mm_conv = .862; %mm per px
nMice = length(perMouseData);
nRows = ceil(sqrt(nMice));
fh = figure; hold on;
dl = [-1, 1];
wind = -30:30;

figure;
for ii = 1:nMice
    % free the data from the structure
    rew_free = perMouseData(ii).turning_traj(ctl_trials{ii});
    dirs = perMouseData(ii).turning_dir(ctl_trials{ii});
    
    free = [];
    for kk = 1:length(rew_free)
        free = cat(2, free, rew_free{kk});
    end
    
    dir_arr = cell2mat(dirs);
    mean_turns = [];
    for jj = 1:length(dl)
        sel = (dir_arr == dl(jj))';
        mean_turns(:,jj) = nanmean(free(:,sel),2);
    end
    subplot(nRows, nMice/nRows, ii);
    plot(mean_turns, 'k'); hold on;
    
    %rew_free2 = perMouseData(ii).turning_traj(ctl2_trials{ii});
    %free2_dir = perMouseData(ii).turning_traj(ctl2_trials{ii});
    
    % Right Occlusions 
    rew_occr = perMouseData(ii).turning_traj(occr_trials{ii});
    dirs = perMouseData(ii).turning_dir(occr_trials{ii});
    
    occr = [];
    for kk = 1:length(rew_occr)
        occr = cat(2, occr, rew_occr{kk});
    end
    
    dir_arr = cell2mat(dirs);
    mean_turns = [];
    for jj = 1:length(dl)
        sel = (dir_arr == dl(jj))';
        mean_turns(:,jj) = nanmean(occr(:,sel),2);
    end
    plot(mean_turns, 'r');
    
    % Left occlusions
    rew_occl = perMouseData(ii).turning_traj(occl_trials{ii});
    dirs = perMouseData(ii).turning_dir(occl_trials{ii});
    
    occl = [];
    for kk = 1:length(rew_occl)
        occl = cat(2, occl, rew_occl{kk});
    end
    
    dir_arr = cell2mat(dirs);
    mean_turns = [];
    for jj = 1:length(dl)
        sel = (dir_arr == dl(jj))';
        mean_turns(:,jj) = nanmean(occl(:,sel),2);
    end
    plot(mean_turns, 'b');
    
%     all_free = cat(1, rew_free, rew_free2);
%     
%     figure(fh);
%     ah = subplot(nRows, nRows, ii); %square, many panels
%     hold on;
%     if(~isempty(rew_occr))
%         plotDistanceHistComparison2(rew_free, rew_occr, following_thresh*mm_conv, '', ah);
%         [p,h] = testCellArrayMedians(all_free, rew_occr);
%         disp(sprintf(' \nWilcoxon test for different medians (right occlusion): h = %f, p = %f\n', h, p)); 
%     end
%     if(~isempty(rew_occl))
%         plotDistanceHistComparison2(rew_free, rew_occl, following_thresh*mm_conv, '', ah);
%         [p,h] = testCellArrayMedians(all_free, rew_occl);
%         disp(sprintf(' \nWilcoxon test for different medians (left occlusion): h = %f, p = %f\n', h, p)); 
%     end
%     if(~isempty(rew_free2))
%         plotDistanceHistComparison2(rew_free, rew_free2, following_thresh*mm_conv, '', ah);
%     end
%     axes(ah); title(mouse_names{ii});
end