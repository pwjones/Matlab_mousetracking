% Process the behavior of the fall cohort of mice.
% Peter Jones, 9.20.2013

fallCohortList; % this script contains the mouse and file names to be included in the analysis.
%spring14CohortList;
base_folder = VIDEO_ROOT;
following_thresh = 20; %mm
clear perMouseData;

[videoList, folder_nums] = listBehavioralVideos(base_folder, folders, mouse_names);

s = matlabpool('size');
if s~=0
    parfor ii = 1:length(videoList)
        perMouseData(ii) = loadMouseTracking(videoList{ii}, 1:length(videoList{ii}), following_thresh);
    end
else
    for ii = 1:length(videoList)
        perMouseData(ii) = loadMouseTracking(videoList{ii}, 1:length(videoList{ii}), following_thresh);
    end
end

%save 'fallCohortData_20mm_followingthresh.mat' perMouseData;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Let's Make a BUNCH OF PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot the time on trail over time.
%plot(1:length(fake_prop), fake_prop *100,'Color', [.25 .25 .25]); hold on;
figure; hold on;
title('Trail Fraction Followed per Trial');
nMice = length(perMouseData);
nRows = ceil(sqrt(nMice));
for ii = 1:nMice
    %subplot(nRows, nRows, ii); %square, many panels
    %nTrials = length(perMouseData(ii).rew_dists);
    nTrails = length(folder_nums{ii}); % added 10/14 - want to plot by day rather than by trial
    
    plot(perMouseData(ii).rew_prop*100, 'g'); hold on;
    plot(perMouseData(ii).dist_prop*100, 'r');
    boxcar = [1 1 1 1 1]./ 5;
    samps_cut = floor(length(boxcar)/2);
    vi = (1+samps_cut):(length(perMouseData(ii).rew_prop)-samps_cut);
    xl = length(vi)+samps_cut;
    %plot([1 xl], 100*[mean(perMouseData(ii).rew_prop) mean(perMouseData(ii).rew_prop)], 'k--');
    %plot([1 xl], [35 35], '--k'); hold on;
    %fake_prop_filt = conv(fake_prop, boxcar,'valid');
    rew_prop_filt = conv(perMouseData(ii).rew_prop, boxcar,'valid');
    dist_prop_filt = conv(perMouseData(ii).dist_prop, boxcar,'valid');
    %vi = (1+samps_cut):(samps_cut+length(rew_prop_filt));
    %flh = plot(vi, fake_prop_filt*100,'LineWidth',2, 'Color', [.25 .25 .25]);
    plot(vi, rew_prop_filt*100, 'g', 'LineWidth',3);
    plot(vi, dist_prop_filt*100, 'r', 'LineWidth',3);
    set(gca, 'TickDir','out', 'fontsize', 16);
    xlabel('Trial #','FontSize', 18);
    ylabel('% Time on Trail','FontSize', 18);
    title(mouse_names{ii});
    %legend({'Rewarded Trail', 'Distracter Trail'});
    %title('Proportion of Time Following Trails');
end

%% Plot the time on trail over time - plot by day
figure; hold on;
title('Trail Fraction Followed per Day');
nMice = length(perMouseData);
nRows = ceil(sqrt(nMice));
max_days = length(folders);
for ii = 1:nMice
    nTrails = length(folder_nums{ii}); % added 10/14 - want to plot by day rather than by trial
    days = sort(unique(folder_nums{ii})); %the day labels
    rew_prop = []; dist_prop = [];
    for jj = 1:length(days)
        sel = folder_nums{ii} == days(jj); 
        rew_prop(jj) = nanmean(perMouseData(ii).rew_prop(sel)) * 100;
        dist_prop(jj) = nanmean(perMouseData(ii).dist_prop(sel)) * 100;
    end
    plot(days, rew_prop, 'g-', 'LineWidth', 1); hold on;
    plot(days, dist_prop, 'r-', 'LineWidth', 1);
    set(gca, 'TickDir','out', 'fontsize', 16);
    xlabel('Training Day','FontSize', 18);
    ylabel('% Time on Trail','FontSize', 18);
    title(mouse_names{ii});
    
    all_days{ii} = days;
    all_rew_prop{ii} = rew_prop;
    all_dist_prop{ii} = dist_prop;
end
% Now we need to add an overall mean line
mean_rew_prop = zeros(max_days, 1)*NaN; mean_dist_prop = zeros(max_days, 1)*NaN;
for ii = 1:max_days
   rew_prop = zeros(nMice,1) * NaN; dist_prop = zeros(nMice,1) * NaN;
   for jj = 1:nMice
      fi = find(all_days{jj} == ii); % find the matching day index for this mouse
      rew_prop(jj) = all_rew_prop{jj}(fi); % select the rewarded prop corresponding to the day
      dist_prop(jj) = all_dist_prop{jj}(fi); 
   end
   mean_rew_prop(ii) = nanmean(rew_prop);
   mean_dist_prop(ii) = nanmean(dist_prop);
end
plot(1:max_days, mean_rew_prop, 'g-', 'LineWidth', 2);
plot(1:max_days, mean_dist_prop, 'r-', 'LineWidth', 2);



%% Let's come up with some distance dependent measures
for jj = 1:length(perMouseData)
    nfiles = length(perMouseData(jj).rew_prop);
    vid_duration = perMouseData(jj).total_frames./ perMouseData(jj).frame_rate;
    perMouseData(jj).max_dist = NaN*zeros(nfiles,2);
    perMouseData(jj).med_dist = NaN*zeros(nfiles,2);
    perMouseData(jj).ncrossings = NaN*zeros(nfiles,2);
    perMouseData(jj).crossing_rates = NaN*zeros(nfiles,2);
    for ii = 1:nfiles
        rdists = perMouseData(jj).rew_dists{ii}; 
        nzi = rdists ~= 0; %these are single frame entries onto the trail
        rdists = rdists(nzi); %eliminate them
        if isempty(rdists) rdists = 0; end %but having an empty variable isn't good later
        ddists = perMouseData(jj).distract_dists{ii};
        nzi = ddists ~= 0;
        ddists = ddists(nzi); 
        if isempty(ddists) ddists = 0; end
        perMouseData(jj).ncrossings(ii,:) = [length(rdists) length(ddists)];
        %rew_crossing_rate(ii) = length(rew_dists{ii})/vid_duration;
        %dist_crossing_rate(ii) = length(distract_dists{ii})/vid_duration;
        perMouseData(jj).max_dist(ii,:) = [max(rdists) max(ddists)];
        perMouseData(jj).med_dist(ii,:) = [median(rdists) median(ddists)];
    end
    perMouseData(jj).crossing_rates(:,1) = perMouseData(jj).ncrossings(:,1)./vid_duration;
    perMouseData(jj).crossing_rates(:,2) = perMouseData(jj).ncrossings(:,2)./vid_duration;
    perMouseData(jj).rew_crossingrates_filt = conv(perMouseData(jj).crossing_rates(:,1), boxcar,'valid');
    perMouseData(jj).dist_crossingrates_filt = conv(perMouseData(jj).crossing_rates(:,2), boxcar,'valid');
    
    % cell arrays of these things
    
end

%% Plotting the distance measures
figure; hold on;
title('Trail Crossings', 'FontSize', 16);
for jj = 1:length(perMouseData)
    hold on;
    plot(perMouseData(jj).crossing_rates(:,1), 'g', 'LineWidth',.5); hold on; 
    plot(perMouseData(jj).crossing_rates(:,2), 'r', 'LineWidth',.5);
    rew_crossingrates_filt = conv(perMouseData(jj).crossing_rates(:,1), boxcar,'valid');
    dist_crossingrates_filt = conv(perMouseData(jj).crossing_rates(:,2), boxcar,'valid');
    vi = (1+samps_cut):(samps_cut+length(rew_crossingrates_filt));
    plot(vi, rew_crossingrates_filt, 'g', 'LineWidth',2); % plot boxcar averaged ones 
    plot(vi, dist_crossingrates_filt, 'r', 'LineWidth',2);
    legend({'Rewarded Trail', 'Distracter Trail'});
    xlabel('Trial #', 'FontSize', 14); ylabel('Trail Crossings per second', 'FontSize', 14);
end

%median distances
figure; hold on;
for jj = 1:length(perMouseData)
    med_dist_filt = [conv( perMouseData(jj).med_dist(:,1), boxcar(:),'valid') conv( perMouseData(jj).med_dist(:,2), boxcar(:),'valid')]; 
    plot(perMouseData(jj).med_dist(:,1), 'g', 'LineWidth', .5); hold on;
    plot(perMouseData(jj).med_dist(:,2), 'r', 'LineWidth', .5);
    vi = (1+samps_cut):(samps_cut+size(med_dist_filt,1));
    plot(vi, med_dist_filt(:,1), 'g', 'LineWidth', 2); 
    plot(vi, med_dist_filt(:,2), 'r', 'LineWidth', 2);    
end
legend({'Rewarded Trail', 'Distracter Trail'});
xlabel('Trial #', 'FontSize', 14); ylabel('Median Following Distance (px)', 'FontSize', 14);

%% Plotting a grid of following distances
figure; hold on;
nMice = length(perMouseData);
nRows = ceil(sqrt(nMice));
for ii=1:nMice
   subplot(nRows, nRows, ii); %square, many panels
   nTrials = length(perMouseData(ii).rew_dists);
   for jj = 1:nTrials
      followingDists = perMouseData(ii).rew_dists{jj};
      x = jj * ones(length(followingDists),1);
      plot(x,followingDists, 'go'); hold on;
      followingDists = perMouseData(ii).distract_dists{jj};
      x = jj * ones(length(followingDists),1);
      plot(x,followingDists, 'rx');
      plot(perMouseData(ii).med_dist(:,1), 'g', 'LineWidth', 2); 
      plot(perMouseData(ii).med_dist(:,2), 'r', 'LineWidth', 2);
   end
   ylim([0 500]);
end

%% Plotting the sum of the distance traveled near the trail
figure; hold on;
title('Total distance traveled near trail');
nMice = length(perMouseData);
nRows = ceil(sqrt(nMice));
for ii=1:4
   %subplot(nRows, nRows, ii); %square, many panels
   nTrials = length(perMouseData(ii).rew_dists);
   rew_dist = []; distract_dist = [];
   for jj = 1:nTrials
      followingDists = perMouseData(ii).rew_dists{jj};
      rew_dist(jj) = sum(followingDists);
      followingDists = perMouseData(ii).distract_dists{jj};
      distract_dist(jj) = sum(followingDists);   
   end
   x = 1:nTrials;
   tt = train_trials{ii};
   boxcar = [1 1 1 1 1]./ 5;
   samps_cut = floor(length(boxcar)/2);
   vi = (1+samps_cut):(length(rew_dist)-samps_cut);
   rew_dist_filt = conv(rew_dist, boxcar,'valid');
   distract_dist_filt = conv(distract_dist, boxcar,'valid');
   plot(x,rew_dist, 'g-', 'LineWidth', .5); hold on;
   plot(x,distract_dist, 'r-', 'LineWidth', .5); hold on;
   plot(vi,rew_dist_filt, 'g-', 'LineWidth', 3); hold on;
   plot(vi,distract_dist_filt, 'r-', 'LineWidth', 3); hold on;
   xlabel('Trial Number'); ylabel('Distance traveled near trail');
   title(mouse_names{ii});
   %ylim([0 1000]);
end


%% Plotting the proportion of the trail traveled per trial
figure; hold on;
title('Trail Fraction Followed per Trial');
nMice = length(perMouseData);
nRows = ceil(sqrt(nMice));
for ii=1:nMice
   %subplot(nRows, nRows, ii); %square, many panels
   nTrials(ii) = length(perMouseData(ii).rew_dists);
   rew_dist = []; distract_dist = [];
   for jj = 1:nTrials(ii)
      %followingDists = perMouseData(ii).rew_propFollowed(jj);
      rew_dist(jj) = perMouseData(ii).rew_propFollowed(jj);
      %followingDists = perMouseData(ii).distract_dists{jj};
      distract_dist(jj) = perMouseData(ii).dist_propFollowed(jj);   
   end
   x = 1:nTrials(ii);
   tt = train_trials{ii};
   plot(x(tt),rew_dist(tt), 'g-', 'LineWidth', .5); hold on;
   plot(x(tt),distract_dist(tt), 'r-', 'LineWidth', .5); hold on;
   boxcar = [1 1 1 1 1]./ 5;
   samps_cut = floor(length(boxcar)/2);
   vi = (1+samps_cut):(length(rew_dist)-samps_cut);
   rew_dist_filt = conv(rew_dist, boxcar,'valid');
   distract_dist_filt = conv(distract_dist, boxcar,'valid');
   plot(vi(tt),rew_dist_filt(tt), 'g-', 'LineWidth', 3); hold on;
   plot(vi(tt),distract_dist_filt(tt), 'r-', 'LineWidth', 3); hold on;
   xlabel('Trial Number','FontSize', 14); ylabel('Prop of Trail Followed','FontSize', 14);
   set(gca, 'TickDir', 'out', 'FontSize', 12)
   title(mouse_names{ii});
   %ylim([0 1000]);
end
maxTrials = max(nTrials);

%%
% Plotting the number of trail pixels followed
figure; hold on;
title('Trail Areas');
nMice = length(perMouseData);
nRows = ceil(sqrt(nMice));
avgRew = NaN.*zeros(maxTrials, nMice);
avgDistract = NaN.*zeros(maxTrials, nMice);
for ii=1:nMice
   %subplot(nRows, nRows, ii); %square, many panels
   nTrials = length(perMouseData(ii).rew_dists);
   rew_dist = []; distract_dist = [];
   for jj = 1:nTrials
      %followingDists = perMouseData(ii).rew_propFollowed(jj);
      rew_dist(jj) = perMouseData(ii).rew_trail_area(jj)*perMouseData(ii).rew_propFollowed(jj);
      %followingDists = perMouseData(ii).distract_dists{jj};
      distract_dist(jj) = perMouseData(ii).dist_trail_area(jj)*perMouseData(ii).dist_propFollowed(jj);   
   end
   x = 1:nTrials;
   tt = train_trials{ii};
   plot(x(tt),rew_dist(tt)/5, 'g-', 'LineWidth', .5); hold on;
   plot(x(tt),distract_dist(tt)/5, 'r-', 'LineWidth', .5); hold on;
   avgRew(1:length(rew_dist),ii) = rew_dist(:);
   avgDistract(1:length(distract_dist),ii) = distract_dist(:);
   vi = (1+samps_cut):(length(rew_dist)-samps_cut);
   rew_dist_filt = conv(rew_dist, boxcar,'valid');
   distract_dist_filt = conv(distract_dist, boxcar,'valid');
   plot(vi(tt),rew_dist_filt(tt)/5, 'g-', 'LineWidth', 3); hold on;
   plot(vi(tt),distract_dist_filt(tt)/5, 'r-', 'LineWidth', 3); hold on;
   xlabel('Trial Number','FontSize', 14); ylabel('Trail Area Followed, mm^2','FontSize', 14);
   set(gca, 'TickDir', 'out', 'FontSize', 12)
   title(mouse_names{ii});
   ylim([0 1500]);
end
%plot(1:maxTrials, nanmean(avgRew,2),'g-', 'LineWidth',2);
%plot(1:maxTrials, nanmean(avgDistract,2),'r-', 'LineWidth',2);


%%
% Plotting the area of trail followed as a rate - This gives a measure of diligence in trail following
% This is the most similar to percentage of time on trail, but is based specifically on the amount of
% trail the animal covers.
figure; hold on;
title('Trail Areas');
nMice = length(perMouseData);
nRows = ceil(sqrt(nMice));
for ii=1:nMice
   %subplot(nRows, nRows, ii); %square, many panels
   nTrials = length(perMouseData(ii).rew_dists);
   rew_dist = []; distract_dist = []; movie_time = [];
   for jj = 1:nTrials
      %followingDists = perMouseData(ii).rew_propFollowed(jj);
      movie_time(jj) = perMouseData(ii).total_frames(jj) ./ perMouseData(ii).frame_rate(jj);
      rew_dist(jj) = perMouseData(ii).rew_trail_area(jj)*perMouseData(ii).rew_propFollowed(jj);
      %followingDists = perMouseData(ii).distract_dists{jj};
      distract_dist(jj) = perMouseData(ii).dist_trail_area(jj)*perMouseData(ii).dist_propFollowed(jj);   
   end
   x = 1:nTrials;
   tt = train_trials{ii};
   rew_dist = rew_dist./movie_time; distract_dist = distract_dist./movie_time; %Make them rates
   plot(x(tt),rew_dist(tt), 'g-', 'LineWidth', .5); hold on;
   plot(x(tt),distract_dist(tt), 'r-', 'LineWidth', .5); hold on;
   vi = (1+samps_cut):(length(rew_dist)-samps_cut);
   rew_dist_filt = conv(rew_dist, boxcar,'valid');
   distract_dist_filt = conv(distract_dist, boxcar,'valid');
   plot(vi(tt),rew_dist_filt(tt), 'g-', 'LineWidth', 3); hold on;
   plot(vi(tt),distract_dist_filt(tt), 'r-', 'LineWidth', 3); hold on;
   xlabel('Trial Number','FontSize', 14); ylabel('Trail Area Following Rate, mm^2/sec','FontSize', 14);
   title(mouse_names{ii});
   set(gca, 'TickDir', 'out', 'FontSize', 12)
   %ylim([0 1000]);
end

%% Plot the following biases of the mice, horizontal position relative to the trail

% We need to define the trials for which the mice are occluded. Since they are occluded in different
% orders this needs to be specified on a per mouse basis.
%ctl_trials = {42:72, 46:76, 43:73, 44:74, 26:45, 26:45, 26:45, 26:45};
%occr_trials = {73:103, 77:107, 74:111, 75:104, 26:35, 26:36, 26:35, 26:35}; %each cell is a mouse
%occl_trials = {16:25, 16:25, 16:25, 16:25, 16:25, 16:25, 16:25, 16:25};
%ctl2_trails = {100:103
nMice = length(perMouseData);
nRows = ceil(sqrt(nMice));
fh = figure; hold on;
for ii = 1:nMice
    rew_free = perMouseData(ii).rew_dists_from_trail_persect(ctl_trials{ii});
    dist_free = perMouseData(ii).distract_dists_from_trail_persect(ctl_trials{ii});
    rew_occ = perMouseData(ii).rew_dists_from_trail_persect(occr_trials{ii});
    dist_occ = perMouseData(ii).distract_dists_from_trail_persect(occr_trials{ii});
    
    figure(fh);
    ah = subplot(nRows, nRows, ii); %square, many panels
    plotDistanceHistComparison(rew_free, dist_free, rew_occ, dist_occ, following_thresh, '', ah);
    axes(ah); title(mouse_names{ii});
end

%% Plot the average nose trajectories

% nMice = length(perMouseData);
% nRows = ceil(sqrt(nMice));
% fh = figure; hold on;
% for ii = 1:nMice
%     traj_window = perMouseData.traj_window;
%     nTrials = length(perMouseData(ii).noseTrajectories);
%     pos_traj = []; neg_traj = [];
%     for jj = 1:nTrials
%         traj = perMouseData(ii).noseTrajectories{jj};
%         pos = perMouseData(ii).traj_dir{jj} == 1;
%         pos_traj = cat(pos_traj
%     end
% end
