% Plot dataset where we've had the mice follow trails of various concentrations

%% Sets vars and loads the dataset %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%concentrationDataList2x; % this script contains the mouse and file names to be included in the analysis.
%concentrationDataList10x;
concentrationDataList3;
%cntnapDataList;
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
%popNoseTraj(perMouseData(1:3), ctl_trials(1:3), {[0 0 0], [1 0 0], [0 0 1]}, []);
%popNoseTraj(perMouseData(4:6), ctl_trials(4:6), {[0 0 0], [1 0 0], [0 0 1]}, []);

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
nMice = length(perMouseData)/nConc;
nRows = ceil(sqrt(nMice));
filtn = 3; boxcar = ones(1,filtn)./filtn; %define the averaging filter kernel
for jj = 1:nMice
    for kk=1:nConc
        ii= (jj-1)*nConc + (kk-1) + 1;
        %subplot(nRows, nRows, ii); %square, many panels
        nTrials = length(perMouseData(ii).rew_dists);
        plot(perMouseData(ii).rew_prop*100, 'g'); hold on;
        plot(perMouseData(ii).dist_prop*100, 'r');
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
end
%% Plot the time on trail over time.
figure; ah1 = axes; hold on;
xlabel('Trial #','FontSize', 18);
ylabel('% Time Following Rewarded Trail','FontSize', 18);
figure; ah2 = axes; hold on;
xlabel('Trial #','FontSize', 18);
ylabel('% Time on Trail','FontSize', 18);
figure; ah3 = axes; hold on;
xlabel('Concentration','FontSize', 18);
ylabel('Prop Following Rewarded','FontSize', 18);
nMice = length(perMouseData)/nConc;
nRows = ceil(sqrt(nMice));
filtn = 3; boxcar = ones(1,filtn)./filtn; %define the averaging filter kernel
mean_prop_follow = []; std_prop_follow = []; n_prop_follow = [];
cl = {[0 0 1], [0 1 0], [0 0 0]}; offset = [-.05 .05];
for jj = 1:nMice
    for kk=1:nConc
        ii= (jj-1)*nConc + (kk-1) + 1;
        %subplot(nRows, nRows, ii); %square, many panels
        nTrials = length(perMouseData(ii).rew_dists);
        plot(ah1, perMouseData(ii).rew_prop*100, 'Color', [0 .8 0]); hold on;
        plot(ah1, perMouseData(ii).dist_prop*100, 'r'); hold on;
        % Let's also do a slightly different measure - ratio of following rewarded/total
        propFollowing = perMouseData(ii).rew_prop./(perMouseData(ii).rew_prop+perMouseData(ii).dist_prop);
        mean_prop_follow(kk,jj) = nanmean(propFollowing); std_prop_follow(kk,jj) = nanstd(propFollowing);
        n_prop_follow(kk,jj) = nTrials;
        plot(ah2, propFollowing*100, 'Color', [0 .8 0]); hold on;
        plot(ah3, (kk+offset(jj))*ones(nTrials,1), propFollowing, '.', 'Color', cl{jj});
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
        plot(ah1, vi, rew_prop_filt*100, 'Color', [0 .8 0], 'LineWidth',3);
        plot(ah1, vi, dist_prop_filt*100, 'r', 'LineWidth',3);
        set(ah1, 'TickDir','out', 'fontsize', 16);
        plot(ah2, vi(:), (rew_prop_filt./(rew_prop_filt+dist_prop_filt))*100, 'Color', 'g', 'LineWidth',3);
        set(ah2, 'TickDir','out', 'fontsize', 16);
    end
end
% Let's produce a plot of the proportion of time on trails spent following the rewarded trail, as a function of concentration

for jj=1:nMice
    plot(ah3, mean_prop_follow(:,jj),'o-', 'Color', cl{jj}); hold on;
    addErrorBars(ah3, 1:nConc, mean_prop_follow(:,jj), std_prop_follow(:,jj)./sqrt(n_prop_follow(:,jj)), cl{jj}, .1);
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
   plot(vi, rew_dist_filt, 'g-', 'LineWidth', 3); hold on;
   plot(vi, distract_dist_filt, 'r-', 'LineWidth', 3); hold on;
   xlabel('Trial Number','FontSize', 14); ylabel('Trail Area Following Rate, mm^2/sec','FontSize', 14);
   title(mouse_names{ii});
   set(gca, 'TickDir', 'out', 'FontSize', 12)
   %ylim([0 1000]);
end

%% Let's calculate some distance dependent measures
for jj = 1:length(perMouseData)
    nfiles = length(perMouseData(jj).rew_prop);
    vid_duration = perMouseData(jj).total_frames./ perMouseData(jj).frame_rate;
    perMouseData(jj).max_dist = NaN*zeros(nfiles,2);
    perMouseData(jj).med_dist = NaN*zeros(nfiles,2);
    perMouseData(jj).ncrossings = NaN*zeros(nfiles,2);
    perMouseData(jj).crossing_rates = NaN*zeros(nfiles,2);
    for ii = 1:nfiles
        rdists = perMouseData(jj).rew_dists{ii}; %Distance traveled near the trail
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

%% %median distances
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

%% Mouse velocities as they are following the trail
mm_conv = 1.16; %mm/px linear
figure; nva = axes; hold on;
%figure; bva = axes; hold on;
colors = {[0 0 0], [1 0 0], [0 0 1], [0 1 0]};
for jj = 1:length(perMouseData)
    bodyVel = []; noseVel = [];
    for kk = 1:length(perMouseData(jj).body_vel)
        for ll= 1:length(perMouseData(jj).body_vel{kk})
            temp = perMouseData(jj).body_vel{kk}{ll};
            bodyVel = cat(1, bodyVel, temp(:));
            temp = perMouseData(jj).nose_vel{kk}{ll};
            noseVel = cat(1, noseVel, temp(:));
        end
    end
    noseVel = noseVel * mm_conv * 60; %conversion: px/frame * mm/px * frames/sec = mm/sec
    velThresh = 10;
    fi = noseVel >= velThresh;
    noseVel = noseVel(fi);
    %[yv, bins] = hist(bodyVel, 200);
    %plot(bva, bins, yv, 'Color', colors{mod(jj-1,length(colors)-1)+1}, 'LineWidth',2);
    %[yv, bins] = hist(noseVel, 200);
    %plot(nva, bins, yv, 'Color', colors{mod(jj-1,length(colors)-1)+1}, 'LineWidth',2);
    %axes(bva);
    %plotEmpiricalCDF({bodyVel}, .1, {colors{mod(jj-1,length(colors)-1)+1}}, '-', bva);
    axes(nva);
    plotEmpiricalCDF({noseVel}, 1, {colors{mod(jj-1,length(colors)-1)+1}}, '-', nva);
    %plot(bva, [1 2 3 4 5; 1 2 3 4 5], [0 0 0 0 0; 1 1 1 1 1], '--k');
end
%plot(bva, [1 2 3 4 5; 1 2 3 4 5], [0 0 0 0 0; 1000 1000 1000 1000 1000], '--k');
%plot(nva, [1 2 3 4 5; 1 2 3 4 5], [0 0 0 0 0; 1000 1000 1000 1000 1000], '--k');