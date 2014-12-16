%% Efficiency by day
%  Plotting the area of trail followed as a rate - This gives a measure of diligence in trail following
%  This is the most similar to percentage of time on trail, but is based specifically on the amount of
%  trail the animal covers. Sort the data by the behavioral day.

% Expect there to be 2 structures present
% 1) a structure caled perMouseData that has a length of the number of mice in the analysis. It contains
% the data
% 2) day_nums - a cell array that gives the days that each trial happens.  Necessary to split trials up
% into days


nMice = length(perMouseData);
nRows = ceil(sqrt(nMice));
max_days = length(folders);
max_trials = 10; %maximum # of trials per day

all_rates = NaN*zeros(nMice, max_trials, max_days, 2); %Dims: 1-Mice, 2-Trials, 3-Days per catagory, 4-Rew/Dist 5-Catagory

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
    daily_rew_dist = NaN*ones(ndays,max_trials); daily_distract_dist = NaN*ones(ndays,max_trials); % dims: days, trials
    
    for jj = 1:length(unique_days)
        day = unique_days(jj);
        sel = find(day == folder_nums{ii});
        ntrial = length(sel);
        all_rates(ii,1:ntrial, jj, 1) = rew_dist(sel);
        all_rates(ii,1:ntrial, jj, 2) = distract_dist(sel);
    end
    
        
end

fh = figure;

mean_individual_mice = squeeze(nanmean(all_rates, 2));
colors = {'g','r'};
for ii = 1:2
    plot(squeeze(mean_individual_mice(:,:,ii))','Color', colors{ii}); hold on;
end

pooled = reshape(all_rates, nMice*max_trials, max_days,2);
mean_pooled = squeeze(nanmean(pooled,1));
for ii = 1:2
    plot(squeeze(mean_pooled(:,ii)),'Color', colors{ii},'LineWidth', 2); hold on;
end

title('Trail Areas');
xlim([1 Inf]);
xlabel('# Days');
ylabel('Avg Trail Following Efficiency (mm/sec)');
ylim([0 60]);