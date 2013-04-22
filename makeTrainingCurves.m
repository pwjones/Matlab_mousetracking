%let's read the file telling us what sections to load

%Analysis Parameters
following_thresh = 15; %px, the distance from the trail the animal can get before it's counted as not following
%vid_list = '3082files.txt';
%vid_list = '9086files.txt';
vid_list = '3090files.txt';
%vid_list = '9085_unbaited_files.txt';
%vid_list = '3083files.txt';
%vid_list = 'video_list.txt';
%hist_trial_range = [1:10; 60:69];
%hist_trial_range = [1:20; 115:134];
hist_trial_range = [2:7; 2:7];
hist_trial_range = [35:49; 20:34];

%%%%%%%%%%%%%%%% Start in on doing things  %%%%%%%%%%%%%%%%
fid = fopen(vid_list, 'r');
res = textscan(fid, '%s [%d %d]');
fnames = res{1};
starts = res{2};
ends = res{3};

base_fname = '/Users/pwjones/Movies/mouse_training/';

% fnames = fnames(4:7);
% starts = starts(4:7);
% ends = ends(4:7);
%fnames = fnames(1:10);
%starts = starts(1:10);
%ends = ends(1:10);
%fnames = fnames(5:end);
%starts = starts(5:end);
%ends = ends(5:end);


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
%    mt.plotFollowing([], following_thresh, 0);
    rew_prop(ii) = mt.propTimeOnTrail([],1,following_thresh);
    dist_prop(ii) = mt.propTimeOnTrail([],2,following_thresh);
    rew_dists{ii} = mt.distanceOnTrail([],1,following_thresh);
    distract_dists{ii} = mt.distanceOnTrail([],2,following_thresh);
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
med_dist = NaN*zeros(nfiles,2);
ncrossings = NaN*zeros(nfiles,2);
crossing_rates = NaN*zeros(nfiles,2);
for ii = 1:nfiles
    rdists = rew_dists{ii}; 
    nzi = rdists ~= 0; %these are single frame entries onto the trail
    rdists = rdists(nzi); %eliminate them
    if isempty(rdists) rdists = 0; end %but having an empty variable isn't good later
    ddists = distract_dists{ii};
    nzi = ddists ~= 0;
    ddists = ddists(nzi); 
    if isempty(ddists) ddists = 0; end
    ncrossings(ii,:) = [length(rdists) length(ddists)];
    %rew_crossing_rate(ii) = length(rew_dists{ii})/vid_duration;
    %dist_crossing_rate(ii) = length(distract_dists{ii})/vid_duration;
    max_dist(ii,:) = [max(rdists) max(ddists)];
    med_dist(ii,:) = [median(rdists) median(ddists)];
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

%median distances
med_dist_filt = [conv(med_dist(:,1), boxcar(:),'valid') conv(med_dist(:,2), boxcar(:),'valid')]; 
figure;
plot(med_dist(:,1), 'g', 'LineWidth', .5); hold on;
plot(med_dist(:,2), 'r', 'LineWidth', .5);
plot(vi, med_dist_filt(:,1), 'g', 'LineWidth', 2); 
plot(vi, med_dist_filt(:,2), 'r', 'LineWidth', 2);
legend({'Rewarded Trail', 'Distracter Trail'});
xlabel('Trial #', 'FontSize', 14); ylabel('Median Following Distance (px)', 'FontSize', 14);

%% 
distance_comp = cell(2,2); %want to do a reward/distractor, early/late comparison
for ii = 1:2
    trials = hist_trial_range(ii,:);
    for jj=1:length(trials)
        rdists = rew_dists{trials(jj)}; 
        nzi = rdists ~= 0; %these are single frame entries onto the trail
        rdists = rdists(nzi); %eliminate them
        
        ddists = distract_dists{trials(jj)};
        nzi = ddists ~= 0;
        ddists = ddists(nzi); 
        
        distance_comp{ii,1} = cat(1, distance_comp{ii,1}, rdists);
        distance_comp{ii,2} = cat(1, distance_comp{ii,2}, ddists);
    end
end
    
counts = cell(2,2); %want to do a reward/distractor, early/late comparison
xbins = 0:5:400;
[counts{1,1}] = hist(distance_comp{1,1}, xbins); counts{1,2} = hist(distance_comp{1,2},xbins);
counts{2,1} = hist(distance_comp{2,1},xbins); counts{2,2} = hist(distance_comp{2,2}, xbins);

% plotting things
xl = [-5 300];  %limits
yl = max([counts{1,1} counts{1,2}]); yl = [-(yl+5) yl+5];
dg = [0 .8 0]; %darker green
dr = [.8 0 0]; %darker red
figure;
% Early Trials
subplot(2,1,1); hold on;
bar(xbins, counts{1,1}, 'FaceColor', dg); hold on;
xmed = median(distance_comp{1,1});
line([xmed xmed], [0 yl(2)], 'Color', dg, 'LineStyle', '--'); 
text(xmed, yl(2)-10, ['median: ' num2str(xmed)], 'color', dg);

bar(xbins, -counts{1,2}, 'FaceColor', dr);
xmed = median(distance_comp{1,2});
line([xmed xmed], [0 yl(1)], 'Color', dr, 'LineStyle', '--'); 
text(xmed, yl(1)+10, ['median: ' num2str(xmed)], 'color', dr);
title('First 8 Trials');
set(gca, 'TickDir', 'out');
xlim(xl); ylim(yl);
% Late Trials
yl = max([counts{2,1} counts{2,2}]);  yl = [-(yl+5) yl+5];
subplot(2,1,2); hold on;
bar(xbins, counts{2,1}, 'FaceColor', dg); hold on;
xmed = median(distance_comp{2,1});
line([xmed xmed], [0 yl(2)], 'Color', dg, 'LineStyle', '--'); 
text(xmed, yl(2)-10, ['median: ' num2str(xmed)], 'color', dg);

bar(xbins, -counts{2,2}, 'FaceColor', dr); hold on;
xmed = median(distance_comp{2,2});
line([xmed xmed], [0 yl(1)], 'Color', dr, 'LineStyle', '--'); 
text(xmed, yl(1)+10, ['median: ' num2str(xmed)], 'color', dr);
xlim(xl); ylim(yl);
title('Last 8 Trials');
xlabel('Following Distance (px)', 'FontSize', 14);
ylabel('# Instances','FontSize', 14);
set(gca, 'TickDir', 'out');
