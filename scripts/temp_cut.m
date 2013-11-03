%% Plot the following biases of the mice, horizontal position relative to the trail

% We need to define the trials for which the mice are occluded. Since they are occluded in different
% orders this needs to be specified on a per mouse basis.
occr_trials = {68:99, 65:77, 26:35, 26:35, 26:35, 26:36, 26:35, 26:35}; %each cell is a mouse
ctl_trials = {39:67, 40:65, 26:45, 26:45, 26:45, 26:45, 26:45, 26:45};
occl_trials = {16:25, 16:25, 16:25, 16:25, 16:25, 16:25, 16:25, 16:25};
%ctl2_trails = {100:103, 
nMice = length(perMouseData);
fh = figure; hold on;
for ii = 1:nMice
    rew_free = perMouseData(ii).rew_dists_from_trail_persect(ctl_trials{ii});
    dist_free = perMouseData(ii).distract_dists_from_trail_persect(ctl_trials{ii});
    rew_occ = perMouseData(ii).rew_dists_from_trail_persect(occr_trials{ii});
    dist_occ = perMouseData(ii).distract_dists_from_trail_persect(occr_trials{ii});
    
    figure(fh);
    ah = subplot(nRows, nRows, ii); %square, many panels
    plotDistanceHistComparison(rew_free, dist_free, rew_occ, dist_occ, 20, '', ah);
    axes(ah); title(mouse_names{ii});
end