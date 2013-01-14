%let's read the file telling us what sections to load
%fid = fopen('short_video_list.txt', 'r');
fid = fopen('video_list.txt', 'r');
res = textscan(fid, '%s [%d %d]');
fnames = res{1};
starts = res{2};
ends = res{3};

base_fname = '/Users/pwjones/Movies/mouse_training/';

% Video properties
%vid = struct('fname', '20121209/9085_2012-12-09-165014-0000.avi', 'timeRange', [5 250]);
%vid(2) = struct('fname', '20121209/9085_2012-12-09-171006-0000.avi', 'timeRange', [8 250]);
nfiles = length(fnames);
rew_prop = NaN*zeros(nfiles,1);
dist_prop = NaN*zeros(nfiles,1);
rew_dists = cell(nfiles,1);
distract_dists = cell(nfiles,1);
total_frames = NaN*zeros(nfiles,1);
frame_rate = NaN*zeros(nfiles,1);

for ii = 1:nfiles
    fullname = fullfile(base_fname, fnames{ii});
    
    mt = MouseTrackerUnder2(fullname, [],[starts(ii) ends(ii)]);
    
    rew_prop(ii) = mt.propTimeOnTrail([],1,5);
    dist_prop(ii) = mt.propTimeOnTrail([],2,5);
    rew_dists{ii} = mt.distanceOnTrail([],1,5);
    distract_dists{ii} = mt.distanceOnTrail([],2,5);
    % now collect a few factors about each of the movies
    total_frames(ii) = mt.nFrames;
    frame_rate(ii) = mt.frameRate;
end

%% Section for plotting curves
figure;
plot(rew_prop, 'g');
hold on;
plot(dist_prop, 'r');
xlabel('Trial #','FontSize', 14);
boxcar = [1 1 1 1 1]./ 5;
samps_cut = floor(length(boxcar)/2);
rew_prop_filt = conv(rew_prop, boxcar,'valid');
dist_prop_filt = conv(dist_prop, boxcar,'valid');
vi = (1+samps_cut):(samps_cut+length(rew_prop_filt));
plot(vi, rew_prop_filt, 'g', 'LineWidth',2);
plot(vi, dist_prop_filt, 'r', 'LineWidth',2);


legend({'Rewarded Trail', 'Distracter Trail'});
title('Proportion of Time Following Trails');

figure; 
plot(rew_prop./(dist_prop+rew_prop), 'k', 'LineWidth',2);
title('Rewarded Following / Total Following');
xlabel('Trial #', 'FontSize', 14);

%% Let's come up with some distance dependent measures
vid_duration = total_frames./ frame_rate;
max_dist = NaN*zeros(nfiles,2);
ncrossings = NaN*zeros(nfiles,2);
crossing_rates = NaN*zeros(nfiles,2);
for ii = 1:nfiles
    ncrossings(ii,:) = [length(rew_dists{ii}) length(distract_dists{ii})];
    %rew_crossing_rate(ii) = length(rew_dists{ii})/vid_duration;
    %dist_crossing_rate(ii) = length(distract_dists{ii})/vid_duration;
    max_dist(ii,:) = [max(rew_dists{ii}) max(distract_dists{ii})];
    med_dist(ii,:) = [median(rew_dists{ii}) median(distract_dists{ii})];
end
crossing_rates(:,1) = ncrossings(:,1)./vid_duration;
crossing_rates(:,2) = ncrossings(:,2)./vid_duration;
rew_crossingrates_filt = conv(crossing_rates(:,1), boxcar,'valid');
dist_crossingrates_filt = conv(crossing_rates(:,2), boxcar,'valid');

figure;
title('Trail Crossings', 'FontSize', 16);
plot(crossing_rates(:,1), 'g', 'LineWidth',.5); hold on; 
plot(crossing_rates(:,2), 'r', 'LineWidth',.5);
plot(vi, rew_crossingrates_filt, 'g', 'LineWidth',2); % plot boxcar averaged ones 
plot(vi, dist_crossingrates_filt, 'r', 'LineWidth',2);
legend({'Rewarded Trail', 'Distracter Trail'});
xlabel('Trial #', 'FontSize', 14); ylabel('Trail Crossings per second', 'FontSize', 14);

figure;
plot(med_dist(:,1), 'g', 'LineWidth', 2); hold on;
plot(med_dist(:,2), 'r', 'LineWidth', 2);
legend({'Rewarded Trail', 'Distracter Trail'});
xlabel('Trial #', 'FontSize', 14); ylabel('Median Following Distance (px)', 'FontSize', 14);

%% 
distance_comp = cell(2,2); %want to do a reward/distractor, early/late comparison
trial_range = [1:8; 37:44];
for ii = 1:2
    trials = trial_range(ii,:);
    for jj=1:length(trials)
        distance_comp{ii,1} = cat(1, distance_comp{ii,1}, rew_dists{trials(jj)});
        distance_comp{ii,2} = cat(1, distance_comp{ii,2}, distract_dists{trials(jj)});
    end
end
    
counts = cell(2,2); %want to do a reward/distractor, early/late comparison
xbins = 0:5:300;
[counts{1,1}] = hist(distance_comp{1,1}, xbins); counts{1,2} = hist(distance_comp{1,2},xbins);
counts{2,1} = hist(distance_comp{2,1},xbins); counts{2,2} = hist(distance_comp{2,2}, xbins);
figure;
subplot(2,1,1); hold on;
bar(xbins, counts{1,1}, 'g'); hold on; 
bar(xbins, -counts{1,2}, 'r'); hold on;
xlim([-5 205]);
ylim([-100 100]);
title('First 8 Trials');
set(gca, 'TickDir', 'out');

subplot(2,1,2); hold on;
bar(xbins, counts{2,1}, 'g'); hold on; 
bar(xbins, -counts{2,2}, 'r'); hold on;
xlim([-5 205]);
ylim([-100 100]);
title('Last 8 Trials');
xlabel('Following Distance (px)', 'FontSize', 14);
ylabel('# Instances','FontSize', 14);
set(gca, 'TickDir', 'out');


   