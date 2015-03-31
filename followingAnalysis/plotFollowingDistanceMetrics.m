function plotFollowingDistanceMetrics(perMouseData, trials, mouse_names)
%
% perMouseData - the data, on a per mouse basis
% trials - the section of the trials data to plot

% Plotting the proportion of the trail traveled per trial
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
   tt = trials{ii};
   x = 1:length(tt);
   plot(x,rew_dist(tt), 'g-', 'LineWidth', .5); hold on;
   plot(x,distract_dist(tt), 'r-', 'LineWidth', .5); hold on;
   %plot([1 length(x)], [mean(rew_dist(tt)), mean(rew_dist(tt))], 'k--', 'LineWidth',1);
   boxcar = [1 1 1 1 1]./ 5;
   samps_cut = floor(length(boxcar)/2);
   vi = (1+samps_cut):(length(rew_dist(tt))-samps_cut);
   rew_dist_filt = conv(rew_dist(tt), boxcar,'valid');
   distract_dist_filt = conv(distract_dist(tt), boxcar,'valid');
   plot(vi,rew_dist_filt, 'g-', 'LineWidth', 3); hold on;
   plot(vi,distract_dist_filt, 'r-', 'LineWidth', 3); hold on;
   xlabel('Trial Number','FontSize', 14); ylabel('Prop of Trail Followed','FontSize', 14);
   set(gca, 'TickDir', 'out', 'FontSize', 12)
   title(mouse_names{ii});
   %ylim([0 1000]);
end
maxTrials = max(nTrials);


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
   tt = trials{ii};
   x = 1:length(tt);
   plot(x,rew_dist(tt), 'g-', 'LineWidth', .5); hold on;
   plot(x,distract_dist(tt), 'r-', 'LineWidth', .5); hold on;
   %plot([1 length(x)], [mean(rew_dist(tt)), mean(rew_dist(tt))], 'k--', 'LineWidth',1);
   %avgRew(1:length(rew_dist),ii) = rew_dist(:);
   %avgDistract(1:length(distract_dist),ii) = distract_dist(:);
   vi = (1+samps_cut):(length(rew_dist(tt))-samps_cut);
   rew_dist_filt = conv(rew_dist(tt), boxcar,'valid');
   distract_dist_filt = conv(distract_dist(tt), boxcar,'valid');
   plot(vi,rew_dist_filt, 'g-', 'LineWidth', 3); hold on;
   plot(vi,distract_dist_filt, 'r-', 'LineWidth', 3); hold on;
   xlabel('Trial Number','FontSize', 14); ylabel('Trail Area Followed, mm^2','FontSize', 14);
   set(gca, 'TickDir', 'out', 'FontSize', 12)
   title(mouse_names{ii});
   ylim([0 10000]);
end
%plot(1:maxTrials, nanmean(avgRew,2),'g-', 'LineWidth',2);
%plot(1:maxTrials, nanmean(avgDistract,2),'r-', 'LineWidth',2);



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
   tt = trials{ii};
   x = 1:length(tt);
   rew_dist = rew_dist./movie_time; distract_dist = distract_dist./movie_time; %Make them rates
   plot(x,rew_dist(tt), 'g-', 'LineWidth', .5); hold on;
   plot(x,distract_dist(tt), 'r-', 'LineWidth', .5); hold on;
   %plot([1 length(x)], [mean(rew_dist(tt)), mean(rew_dist(tt))], 'k--', 'LineWidth',1);
   %avgRew(1:length(rew_dist),ii) = rew_dist(:);
   %avgDistract(1:length(distract_dist),ii) = distract_dist(:);
   vi = (1+samps_cut):(length(rew_dist(tt))-samps_cut);
   rew_dist_filt = conv(rew_dist(tt), boxcar,'valid');
   distract_dist_filt = conv(distract_dist(tt), boxcar,'valid');
   plot(vi,rew_dist_filt, 'g-', 'LineWidth', 3); hold on;
   plot(vi,distract_dist_filt, 'r-', 'LineWidth', 3); hold on;
   xlabel('Trial Number','FontSize', 14); ylabel('Trail Area Following Rate, mm^2/sec','FontSize', 14);
   title(mouse_names{ii});
   set(gca, 'TickDir', 'out', 'FontSize', 12)
   %ylim([0 1000]);
end