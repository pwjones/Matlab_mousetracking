%% Plot the following biases of the mice, horizontal position relative to the trail

% We need to define the trials for which the mice are occluded. Since they are occluded in different
% orders this needs to be specified on a per mouse basis.
%ctl_trials = {42:72, 46:76, 43:73, 44:74, 26:45, 26:45, 26:45, 26:45};
%occr_trials = {73:103, 77:107, 74:111, 75:104, 26:35, 26:36, 26:35, 26:35}; %each cell is a mouse
%occl_trials = {16:25, 16:25, 16:25, 16:25, 16:25, 16:25, 16:25, 16:25};
%ctl2_trails = {100:103
mm_conv = .862; %mm per px
nMice = length(perMouseData);
nRows = ceil(sqrt(nMice));
fh = figure; hold on;
for ii = 1:nMice
    rew_free = perMouseData(ii).rew_dists_from_trail_persect(ctl_trials{ii});
    for jj=1:length(rew_free)
        rew_free{jj} = rew_free{jj};
    end
    rew_free2 = perMouseData(ii).rew_dists_from_trail_persect(ctl2_trials{ii});
    for jj=1:length(rew_free2)
        rew_free2{jj} = rew_free2{jj};
    end
    rew_occr = perMouseData(ii).rew_dists_from_trail_persect(occr_trials{ii});
    for jj=1:length(rew_occr)
        rew_occr{jj} =  rew_occr{jj};
    end
    rew_occl = perMouseData(ii).rew_dists_from_trail_persect(occl_trials{ii});
    for jj=1:length(rew_occl)
        rew_occl{jj} = rew_occl{jj};
    end
    
    all_free = cat(1, rew_free, rew_free2);
    
    figure(fh);
    ah = subplot(nRows, nRows, ii); %square, many panels
    hold on;
    if(~isempty(rew_occr))
        plotDistanceHistComparison2(rew_free, rew_occr, following_thresh*mm_conv, '', ah);
        [p,h] = testCellArrayMedians(all_free, rew_occr);
        disp(sprintf(' \nWilcoxon test for different medians (right occlusion): h = %f, p = %f\n', h, p)); 
    end
    if(~isempty(rew_occl))
        plotDistanceHistComparison2(rew_free, rew_occl, following_thresh*mm_conv, '', ah);
        [p,h] = testCellArrayMedians(all_free, rew_occl);
        disp(sprintf(' \nWilcoxon test for different medians (left occlusion): h = %f, p = %f\n', h, p)); 
    end
    if(~isempty(rew_free2))
        plotDistanceHistComparison2(rew_free, rew_free2, following_thresh*mm_conv, '', ah);
    end
    axes(ah); title(mouse_names{ii});
end