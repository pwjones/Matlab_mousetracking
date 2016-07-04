%%
% Plotting the area of trail followed as a rate - This gives a measure of diligence in trail following
% This is the most similar to percentage of time on trail, but is based specifically on the amount of
% trail the animal covers.
cats = {'train', 'ctl', 'occ1', 'ctl2', 'occ2'};

fh = figure;
title('Trail Areas');
nMice = length(perMouseData);
nRows = ceil(sqrt(nMice));
max_days = length(folders);
max_trials = 10;
all_rates = NaN*zeros(nMice, max_trials, 15, 2, length(cats)); %Dims: 1-Mice, 2-Trials, 3-Days per catagory, 4-Rew/Dist 5-Catagory
for ii=1:nMice
    nTrials = length(perMouseData(ii).rew_dists);
    rew_dist = []; distract_dist = []; movie_time = [];
    abs_days = folder_nums{ii};
    for jj = 1:nTrials
        movie_time(jj) = perMouseData(ii).total_frames(jj) ./ perMouseData(ii).frame_rate(jj);
        rew_dist(jj) = perMouseData(ii).rew_trail_area(jj)*perMouseData(ii).rew_propFollowed(jj);
        distract_dist(jj) = perMouseData(ii).dist_trail_area(jj)*perMouseData(ii).dist_propFollowed(jj);
    end
    rew_dist = rew_dist./movie_time; distract_dist = distract_dist./movie_time; %Make them rates
    
    unique_days = sort(unique(folder_nums{ii}));
    ndays = length(unique_days);
    types = cell(ndays,1);
    daily_rew_dist = NaN*ones(ndays,max_trials); daily_distract_dist = NaN*ones(ndays,max_trials); % dims: days, trials
    
    for jj = 1:length(unique_days)
        day = unique_days(jj);
        sel = find(unique_days(jj) == folder_nums{ii});
        ntrial = length(sel);
        daily_rew_dist(jj, 1:ntrial) = rew_dist(sel);
        daily_distract_dist(jj, 1:ntrial) = distract_dist(sel);
        
        if ~isempty(find(day == train_dirs{ii},1))
            types{jj} = 'train';
        elseif ~isempty(find(day == ctl_dirs{ii},1))
            types{jj} = 'ctl';
        elseif ~isempty(find(day == ctl2_dirs{ii},1))
            types{jj} = 'ctl2';
        elseif ~isempty(find(day == occr_dirs{ii},1))
            types{jj} = 'occ1';
        elseif ~isempty(find(day == occl_dirs{ii},1))
            types{jj} = 'occ2';
        else
            types{jj} = 'other';
        end
        
    end
    
    for cc = 1:length(cats)
        subplot(1, length(cats), cc); hold on;
        sel = strcmp(cats{cc}, types);
        if (sum(sel) > 0)
            plot(nanmean(daily_distract_dist(sel,:),2), 'r');
            plot(nanmean(daily_rew_dist(sel,:),2), 'g');
            nsel = sum(sel);
            % Assign to a multi-dim array to make averages easier
            all_rates(ii,:,1:nsel, 1, cc) = daily_rew_dist(sel,:)';
            all_rates(ii,:,1:nsel, 2, cc) = daily_distract_dist(sel,:)';
        end
        xlim([1 Inf]);
        xlabel('# Days');
        ylabel('Avg Trail Following Efficiency (mm/sec)');
        ylim([0 60]);
        set(gca, 'TickDir', 'out');  
    end
end
% Assemble and plot the mean rates
arsiz = size(all_rates);
trial_rates = reshape(all_rates, [arsiz(1)*arsiz(2), arsiz(3), arsiz(4:end)]); %
mean_rates = squeeze(nanmean(trial_rates, 1));
for cc = 1:length(cats)
   subplot(1, length(cats), cc);
   plot(mean_rates(:,2,cc), 'r', 'LineWidth', 2);
   plot(mean_rates(:,1,cc), 'g', 'LineWidth', 2); 
end

trsiz = size(trial_rates);
bar_data = reshape(permute(trial_rates, [1 2 4 3]), [trsiz(1)*trsiz(2), trsiz(4), trsiz(3)]);
ntrials = squeeze(sum(isnan(bar_data),1));
bary = squeeze(nanmean(bar_data,1));
bary = bary(2:5,:)';
barx = [.8 1.2 1.8 2.2 2.8 3.2 3.8 4.2];
barerr = squeeze(nanstd(bar_data, 0,1))./ sqrt(ntrials);
barerr = barerr(2:5,:)';
figure; bh = bar(barx, bary(:),'BarWidth', 1);
addErrorBars(gca, barx, bary(:), barerr(:), 'k', .1)
ylabel('Mean Following Efficiency (mm/sec)');
set(gca, 'FontSize', 22)

ctl = squeeze(bar_data(:, 2, 1)); ctl = ctl(~isnan(ctl));
occ = squeeze(bar_data(:, 3, 1)); occ = occ(~isnan(occ));
[p,h] = ranksum(ctl, occ);
disp(sprintf('Wilcoxon, occlusion 1: h = %i, p = %f', h, p));

occ = squeeze(bar_data(:, 5, 1)); occ = occ(~isnan(occ));
[p,h] = ranksum(ctl, occ);
disp(sprintf('Wilcoxon, occlusion 2: h = %i, p = %f', h, p));

ctl2 = squeeze(bar_data(:, 4, 1)); ctl2 = ctl2(~isnan(ctl2));
[p,h] = ranksum(ctl, ctl2);
disp(sprintf('Wilcoxon, control2: h = %i, p = %f', h, p));

figure(fh)
for ii = 1:size(bary,2)
   subplot(1, length(cats), ii+1);
   plot([1,15], [bary(2,ii) bary(2,ii)], 'r--', 'LineWidth', 2);
   plot([1,15], [bary(1,ii), bary(1,ii)], 'g--', 'LineWidth', 2);
    
end
%% Plot the time on trail over time - plot by day
% figure; hold on;
% title('Trail Fraction Followed per Day');
% nMice = length(perMouseData);
% nRows = ceil(sqrt(nMice));
% max_days = length(folders);
% for ii = 1:nMice
%     nTrails = length(folder_nums{ii}); % added 10/14 - want to plot by day rather than by trial
%     days = sort(unique(folder_nums{ii})); %the day ids
%     day_labels = 1:length(days);
%     rew_prop = []; dist_prop = [];
%     for jj = 1:length(days)
%         sel = folder_nums{ii} == days(jj); 
%         rew_prop(jj) = nanmean(perMouseData(ii).rew_prop(sel)) * 100;
%         dist_prop(jj) = nanmean(perMouseData(ii).dist_prop(sel)) * 100;
%     end
%     plot(day_labels, rew_prop, 'g-', 'LineWidth', 1); hold on;
%     plot(day_labels, dist_prop, 'r-', 'LineWidth', 1);
%     set(gca, 'TickDir','out', 'fontsize', 16);
%     xlabel('Training Day','FontSize', 18);
%     ylabel('% Time on Trail','FontSize', 18);
%     title(mouse_names{ii});
%     
%     all_days{ii} = days;
%     all_rew_prop{ii} = rew_prop;
%     all_dist_prop{ii} = dist_prop;
% end
% % Now we need to add an overall mean line
% rew_prop = zeros(nMice,max_days) * NaN; dist_prop = zeros(nMice,max_days) * NaN;
% for ii = 1:nMice
%    n_days = length(all_days{ii});
%    rew_prop(ii, 1:n_days) = all_rew_prop{ii};
%    dist_prop(ii, 1:n_days) = all_dist_prop{ii};
% end
% mean_rew_prop = nanmean(rew_prop,1);
% mean_dist_prop = nanmean(dist_prop,1);
% plot(1:max_days, mean_rew_prop, 'g-', 'LineWidth', 2);
% plot(1:max_days, mean_dist_prop, 'r-', 'LineWidth', 2);