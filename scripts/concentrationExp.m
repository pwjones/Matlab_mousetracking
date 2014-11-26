% Plot dataset where we've had the mice follow trails of various concentrations

%% Sets vars and loads the dataset %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%concentrationDataList2x; % this script contains the mouse and file names to be included in the analysis.
%concentrationDataList10x;
%concentrationDataList3;
concentrationDataList4;
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
%nMice = 2;
%nConc = 3;
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
    plotDistanceHistComparison(rew{1},rew{2}, rew{3}, rew{4}, following_thresh, '', ah);
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
ylabel('% Time on Trail','FontSize', 18);
figure; ah2 = axes; hold on;
xlabel('Trial #','FontSize', 18);
ylabel('% Time Following Rewarded Trail','FontSize', 18);
figure; ah3 = axes; hold on;
xlabel('Concentration','FontSize', 18);
ylabel('Prop Following Rewarded','FontSize', 18);
nMice = length(perMouseData)/nConc;
nRows = ceil(sqrt(nMice));
filtn = 3; boxcar = ones(1,filtn)./filtn; %define the averaging filter kernel
mean_prop_follow = []; std_prop_follow = []; n_prop_follow = [];
cl = {[0 0 1], [0 1 0], [0 0 0]}; offset = [-.05 .05];
rew_colors = blackGradColormap([0 1 0], nConc+1); rew_colors = rew_colors(2:end,:);
dist_colors = blackGradColormap([1 0 0], nConc+1); dist_colors = dist_colors(2:end,:);
for jj = 1:nMice
    for kk=1:nConc
        ii= (jj-1)*nConc + (kk-1) + 1;
        %subplot(nRows, nRows, ii); %square, many panels
        nTrials = length(perMouseData(ii).rew_dists);
        plot(ah1, perMouseData(ii).rew_prop*100, 'Color', rew_colors(kk,:)); hold on;
        plot(ah1, perMouseData(ii).dist_prop*100, 'Color', dist_colors(kk,:)); hold on;
        % Let's also do a slightly different measure - ratio of following rewarded/total
        propTimeFollowingRew = perMouseData(ii).rew_prop./(perMouseData(ii).rew_prop+perMouseData(ii).dist_prop);
        mean_prop_follow(kk,jj) = nanmean(propTimeFollowingRew); std_prop_follow(kk,jj) = nanstd(propTimeFollowingRew);
        n_prop_follow(kk,jj) = nTrials;
        plot(ah2, propTimeFollowingRew*100, 'Color', rew_colors(kk,:)); hold on;
        plot(ah3, (kk+offset(jj))*ones(nTrials,1), propTimeFollowingRew, '.', 'Color', cl{jj});
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
        plot(ah1, vi, rew_prop_filt*100, 'Color', rew_colors(kk,:), 'LineWidth',3);
        plot(ah1, vi, dist_prop_filt*100, 'Color', dist_colors(kk,:), 'LineWidth',3);
        set(ah1, 'TickDir','out', 'fontsize', 16);
        plot(ah2, vi(:), (rew_prop_filt./(rew_prop_filt+dist_prop_filt))*100, 'Color', rew_colors(kk,:), 'LineWidth',3);
        set(ah2, 'TickDir','out', 'fontsize', 16);
    end
end
% Let's produce a plot of the proportion of time on trails spent following the rewarded trail, as a function of concentration

for jj=1:nMice
    plot(ah3, mean_prop_follow(:,jj),'o-', 'Color', cl{jj}); hold on;
    addErrorBars(ah3, 1:nConc, mean_prop_follow(:,jj), std_prop_follow(:,jj)./sqrt(n_prop_follow(:,jj)), cl{jj}, .1);
end

%% Plotting the percent of the trail explored in each case 
maxtrial = 12;
nMice = 2; nConc = 4;
propAreaFollowedRatio = NaN(maxtrial, nMice*nConc);
propTimeFollowedRatio = NaN(maxtrial, nMice*nConc);
for jj = 1:nMice
    for kk=1:nConc
       ii= (jj-1)*nConc + (kk-1) + 1;
       pmd = perMouseData(ii);
       g(ii) = ii; %group variable for box plot
       conc(ii) = kk;
       n = length(pmd.rew_prop);
       propTimeFollowedRatio(1:length(pmd.rew_prop),ii) = pmd.rew_prop ./ (pmd.dist_prop + pmd.rew_prop);
       propAreaFollowedRatio(1:length(pmd.rew_propFollowed),ii) = pmd.rew_propFollowed ./ (pmd.dist_propFollowed + pmd.rew_propFollowed);
    end
end
% plotting measures of the area of trail followed in each case
figure;
boxplot(propAreaFollowedRatio, g, 'notch', 'on');
figure;
conc = reshape(conc, nConc, []);
meanAreaFollowedRatio = reshape(nanmean(propAreaFollowedRatio)', nConc, []);
semAreaFollowedRatio = reshape(nanstd(propAreaFollowedRatio)'./n, nConc, []);
errorbar(conc, meanAreaFollowedRatio, semAreaFollowedRatio);
xlabel('Concentration'); ylabel('Area of Rewared Trial Followed/ Area all trails followed');
%plotting measures of the time following the correct trail
figure;
boxplot(propTimeFollowedRatio, g, 'notch', 'on');
figure;
conc = reshape(conc, nConc, []);
meanTimeFollowedRatio = reshape(nanmean(propTimeFollowedRatio)', nConc, []);
semTimeFollowedRatio = reshape(nanstd(propTimeFollowedRatio)'./n, nConc, []);
errorbar(conc, meanTimeFollowedRatio, semTimeFollowedRatio);
 xlabel('Concentration'); ylabel('% Following time on Rewarded trail');

%% Plotting the area of trail followed as a rate
% This gives a measure of diligence in trail following
% This is the most similar to percentage of time on trail, but is based specifically on the amount of
% trail the animal covers.
figure; hold on;
title('Trail Areas');
rew_colors = blackGradColormap([0 1 0], nConc+1); rew_colors = rew_colors(2:end,:);
dist_colors = blackGradColormap([1 0 0], nConc+1); dist_colors = dist_colors(2:end,:);
for ll = 1:nMice
    for kk=1:nConc
        ii= (ll-1)*nConc + (kk-1) + 1;
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
        plot(x(tt),rew_dist(tt), '-', 'LineWidth', .5, 'Color', rew_colors(kk,:)); hold on;
        plot(x(tt),distract_dist(tt), '-', 'LineWidth', .5, 'Color', dist_colors(kk,:)); hold on;
        vi = (1+samps_cut):(length(rew_dist)-samps_cut);
        rew_dist_filt = conv(rew_dist, boxcar,'valid');
        distract_dist_filt = conv(distract_dist, boxcar,'valid');
        plot(vi, rew_dist_filt, '-', 'LineWidth', 3, 'Color', rew_colors(kk,:)); hold on;
        plot(vi, distract_dist_filt, '-', 'LineWidth', 3, 'Color', dist_colors(kk,:)); hold on;
        xlabel('Trial Number','FontSize', 14); ylabel('Trail Area Following Rate, mm^2/sec','FontSize', 14);
        title(mouse_names{ii});
        set(gca, 'TickDir', 'out', 'FontSize', 12)
        %ylim([0 1000]);
    end
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

%% Mouse velocities as they are following the trail
mm_conv = .862; %mm/px linear
figure; nva = axes; hold on;
max_l = 0;
colors = blackGradColormap([0 1 0], nConc+1); colors = colors(2:end,:);
for jj = 1:length(perMouseData)
    conc_i = mod(jj-1, nConc)+1;
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
    axes(nva);
    plotEmpiricalCDF({noseVel}, 1, {colors(conc_i,:)}, '-', nva);
    allNoseVel{jj} = noseVel;
    max_l = max(length(noseVel), max_l);
    median_noseVel(jj) = nanmedian(noseVel);
    conc(jj) = conc_i;
end
size(median_noseVel)
median_noseVel = reshape(median_noseVel, nConc,[]);
conc = reshape(conc, nConc, []);
figure;
plot(conc, median_noseVel, 'o-');
