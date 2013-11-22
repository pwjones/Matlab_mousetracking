% Make a graph of the median nose position for the mice loaded.

% We need to define the trials for which the mice are occluded. Since they are occluded in different
% orders this needs to be specified on a per mouse basis.

% perMouseData is present and has multiple mice
following_thresh = 20;
mm_conv = 1.16; %mm/px conversion
mouse_names = {'16220', '16221', '16227', '16231', '3090', '3091'};
%mouse_names = {'16220', '16221', '16227', '16231'};
train_trials = {1:72, 1:76, 1:73, 1:74, 1:35, 1:70};
ctl_trials = {42:72, 46:76, 43:73, 44:74, 35:67, 70:100};
occr_trials = {[], 77:107, 74:111, 75:104, 16:34, 101:138}; %each cell is a mouse
occl_trials = {73:103, [], [], [], 68:94, 193:212};
ctl2_trials = {[], [], [], [], [], []};
nMice = length(perMouseData);

fh = figure; hold on;
% want to compile the median nose distance from trail while following
med_nose_pos_ctl = NaN*zeros(nMice, 1);
med_nose_pos_ctl2 = NaN*zeros(nMice, 1);
med_nose_pos_occr = NaN*zeros(nMice, 1);
med_nose_pos_occl = NaN*zeros(nMice, 1);
for ii = 1:nMice
    all_dists = perMouseData(ii).rew_dists_from_trail_persect(ctl_trials{ii});
    following_dists = [];
    for jj = 1:length(all_dists)
        trial_dists = all_dists{jj};
        %fi = find(abs(trial_dists) <= following_thresh & ~isnan(trial_dists));
        fi = find(~isnan(trial_dists));
        following_dists = cat(1, following_dists, trial_dists(fi));
    end
    nose_pos_ctl{ii} = following_dists;
    med_nose_pos_ctl(ii) = median(following_dists);
    
    
    all_dists = perMouseData(ii).rew_dists_from_trail_persect(ctl2_trials{ii});
    following_dists = [];
    for jj = 1:length(all_dists)
        trial_dists = all_dists{jj};
        %fi = find(abs(trial_dists) <= following_thresh & ~isnan(trial_dists));
        fi = find(~isnan(trial_dists));
        following_dists = cat(1, following_dists, trial_dists(fi));
    end
    nose_pos_ctl2{ii} = following_dists;
    med_nose_pos_ctl2(ii) = median(following_dists);
    
    all_dists = perMouseData(ii).rew_dists_from_trail_persect(occr_trials{ii});
    following_dists = [];
    for jj = 1:length(all_dists)
        trial_dists = all_dists{jj};
        %fi = find(abs(trial_dists) <= following_thresh & ~isnan(trial_dists));
        fi = find(~isnan(trial_dists));
        following_dists = cat(1, following_dists, trial_dists(fi));
    end
    nose_pos_occr{ii} = following_dists;
    med_nose_pos_occr(ii) = median(following_dists);
    
    all_dists = perMouseData(ii).rew_dists_from_trail_persect(occl_trials{ii});
    following_dists = [];
    for jj = 1:length(all_dists)
        trial_dists = all_dists{jj};
        %fi = find(abs(trial_dists) <= following_thresh & ~isnan(trial_dists));
        fi = find(~isnan(trial_dists));
        following_dists = cat(1, following_dists, trial_dists(fi));
    end
    nose_pos_occl{ii} = following_dists;
    med_nose_pos_occl(ii) = median(following_dists);
    
    % Let's do some statistical tests on the resulting distributions
    hl = 0; hr = 0;
    if ~isempty(nose_pos_ctl{ii}) && ~isempty(nose_pos_occr{ii})
        [p,hr] = ranksum(nose_pos_ctl{ii}, nose_pos_occr{ii});
        if hr sig = ''; else sig = 'NOT'; end
        sig_str = sprintf('Mouse %s Wilcoxon: Unoccluded-Right Occlusion is %s significantly different: h=%i p=%f', ...
                           mouse_names{ii}, sig, hr, p);
        disp(sig_str);
    end
    
    if ~isempty(nose_pos_ctl{ii}) && ~isempty(nose_pos_occl{ii})
        [p,hl] = ranksum(nose_pos_ctl{ii}, nose_pos_occl{ii});
        if hl sig = ''; else sig = 'NOT'; end
        sig_str = sprintf('Mouse %s Wilcoxon: Unoccluded-Left Occlusion is %s significantly different: h=%i p=%f', ...
                           mouse_names{ii}, sig, hl, p);
        disp(sig_str);
    end
    
    % now let's plot the medians with fill indicating significance
    y = [med_nose_pos_occl(ii), med_nose_pos_ctl(ii), med_nose_pos_occr(ii)]; 
    x = 1:3;
    plot(x, y, 'o-k', 'MarkerSize', 10);
    % Replot the filled markers over if they are significant
    if hl
        plot(1,y(1), 'ok', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
    end
    if hr
        plot(3,y(3), 'ok', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
    end
end

% Make some asthetic changes to the plot
set(gca, 'FontSize', 12, 'xtick', 1:3, 'xticklabel', {'Left Occluded', 'Unoccluded', 'Right Occluded'}, 'TickDir', 'out');
ylabel('Median Distance from Trail (mm)', 'FontSize', 14);

