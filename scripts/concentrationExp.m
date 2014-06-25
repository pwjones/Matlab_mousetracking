% Plot dataset where we've had the mice follow trails of various concentrations

%% Sets vars and loads the dataset %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%concentrationDataList2x; % this script contains the mouse and file names to be included in the analysis.
concentrationDataList10x;
%spring14CohortList;
base_folder = VIDEO_ROOT;
following_thresh = 20; %mm
clear perMouseData;

videoList = listBehavioralVideos(base_folder, folders, mouse_names);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Run analysis on the dataset %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Averaged Trail Crossing Triggered Nose Trajectories
popNoseTraj(perMouseData(1:3), ctl_trials(1:3), {[0 0 0], [1 0 0], [0 0 1]}, []);
popNoseTraj(perMouseData(4:6), ctl_trials(4:6), {[0 0 0], [1 0 0], [0 0 1]}, []);

%% Nose position CDF while following
nMice = 2;
nConc = 3;
nRows = ceil(sqrt(nMice));
fh = figure; hold on;
for jj = 1:nMice
    for kk = 1:nConc
        ii= (jj-1)*nConc + (kk-1) + 1;
        %rew{kk} = perMouseData(ii).rew_dists_from_trail_persect(ctl_trials{ii});
        %dist{kk} = perMouseData(ii).distract_dists_from_trail_persect(ctl_trials{ii}); 
        rew{kk} = perMouseData(ii).rew_dists_from_trail(ctl_trials{ii});
        dist{kk} = perMouseData(ii).distract_dists_from_trail(ctl_trials{ii});
    end
    figure(fh);
    ah = subplot(nRows, nRows, jj); %square, many panels
    plotDistanceHistComparison(rew{1},rew{2}, rew{3}, dist{1}, following_thresh, '', ah);
    axes(ah); title(mouse_names{ii});
end

%% Plot the time on trail over time.
%plot(1:length(fake_prop), fake_prop *100,'Color', [.25 .25 .25]); hold on;
figure; hold on;
title('Trail Fraction Followed per Trial');
nMice = length(perMouseData);
nRows = ceil(sqrt(nMice));
for ii = 1:(nMice*nConc)
    %subplot(nRows, nRows, ii); %square, many panels
    nTrials = length(perMouseData(ii).rew_dists);
    plot(perMouseData(ii).rew_prop*100, 'g'); hold on;
    plot(perMouseData(ii).dist_prop*100, 'r');
    boxcar = [1 1 1]./ 5;
    samps_cut = floor(length(boxcar)/2);
    vi = (1+samps_cut):(length(perMouseData(ii).rew_prop)-samps_cut)
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

%% Plotting the area of trail followed as a rate
% This gives a measure of diligence in trail following
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