%% Plot the following biases of the mice, horizontal position relative to the trail

% We need to define the trials for which the mice are occluded. Since they are occluded in different
% orders this needs to be specified on a per mouse basis.
%ctl_trials = {42:72, 46:76, 43:73, 44:74, 26:45, 26:45, 26:45, 26:45};
%occr_trials = {73:103, 77:107, 74:111, 75:104, 26:35, 26:36, 26:35, 26:35}; %each cell is a mouse
%occl_trials = {16:25, 16:25, 16:25, 16:25, 16:25, 16:25, 16:25, 16:25};
%ctl2_trails = {100:103
mm_conv = 1.16; %mm per px
nMice = length(perMouseData);
nRows = ceil(sqrt(nMice));
fh = figure; hold on;
for ii = 1:nMice
    rew_free = perMouseData(ii).rew_dists_from_trail_persect(ctl_trials{ii});
    for jj=1:length(rew_free)
        rew_free{jj} = mm_conv * rew_free{jj};
    end
    rew_occr = perMouseData(ii).rew_dists_from_trail_persect(occr_trials{ii});
    for jj=1:length(rew_occr)
        rew_occr{jj} = mm_conv * rew_occr{jj};
    end
    rew_occl = perMouseData(ii).rew_dists_from_trail_persect(occl_trials{ii});
    for jj=1:length(rew_occl)
        rew_occl{jj} = mm_conv * rew_occl{jj};
    end
    
    
    figure(fh);
    ah = subplot(nRows, nRows, ii); %square, many panels
    hold on;
    plotDistanceHistComparison2(rew_free, rew_occr, following_thresh, '', ah);
    plotDistanceHistComparison2(rew_free, rew_occl, following_thresh, '', ah);
    axes(ah); title(mouse_names{ii});
end