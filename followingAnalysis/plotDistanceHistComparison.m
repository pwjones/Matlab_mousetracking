function plotDistanceHistComparison(rew_cell, dist_cell, rew_cell2, dist_cell2, dist_thresh, plotOptions, varargin)
% function plotDistanceHists(dist_cell)
%
% Function will plot histograms of the frame by frame distance from the
% path for 
min_dist = 0; %.75;

if ~isempty(varargin)
    cdf_ah = varargin{1};
else
    figure; cdf_ah = axes; hold on;
end
if isempty(plotOptions)
   plotOptions = {'-','-','--','--'}; 
end

nfiles = [length(rew_cell), length(dist_cell); length(rew_cell2), length(dist_cell2)];
rew = {rew_cell, rew_cell2};
dist = {dist_cell, dist_cell2};
distance_comp = cell(2,2); %want to do a reward/distractor, early/late comparison
for ii = 1:2
    %rewarded path
    trials = nfiles(ii,1);
    for jj=1:trials
        %compile distances from the trail for two conditions
        r_trials = rew{ii};
        rdists = r_trials{jj};
        inci = (abs(rdists) <= dist_thresh) & ~isnan(rdists) & (abs(rdists) >= min_dist);
        rdists = rdists(inci);
        distance_comp{ii,1} = cat(1, distance_comp{ii,1}, rdists(:));
    end
    % distractor path
    trials = nfiles(ii,2);
    for jj=1:trials
        d_trials = dist{ii};
        ddists = d_trials{jj};
        inci = (abs(ddists) <= dist_thresh) & ~isnan(ddists) & (abs(ddists) >= min_dist);
        ddists = ddists(inci);
        distance_comp{ii,2} = cat(1, distance_comp{ii,2}, ddists(:));
    end
end
    
counts = cell(2,2); %want to do a reward/distractor, early/late comparison
xbins = -dist_thresh:.5:dist_thresh;
[counts{1,1}] = hist(distance_comp{1,1}, xbins); counts{1,2} = hist(distance_comp{1,2},xbins);
counts{2,1} = hist(distance_comp{2,1},xbins); counts{2,2} = hist(distance_comp{2,2}, xbins);

% plotting things
xl = [-dist_thresh dist_thresh];  %limits
yl = max([counts{1,1} counts{2,1}]); yl = [-(yl+5) yl+5];
dg = [0 .8 0]; %darker green
dr = [.8 0 0]; %darker red
figure;
% Rewarded Path Following
subplot(2,1,1); hold on;
bar(xbins, counts{1,1}, 'FaceColor', dg); hold on;
xmed = double(median(distance_comp{1,1}));
line([xmed xmed], [0 yl(2)], 'Color', dg, 'LineStyle', '--'); 
text(xmed, yl(2), ['median: ' num2str(xmed)], 'Color', dg);

bar(xbins, -counts{2,1}, 'FaceColor', dg); hold on;
xmed = double(median(distance_comp{2,1}));
line([xmed xmed], [0 -yl(2)], 'Color', dg, 'LineStyle', '--'); 
text(xmed, -yl(2)-10, ['median: ' num2str(xmed)], 'color', dg);

% KS test the distributions
[h,p,ks2stat] = kstest2(distance_comp{1,1}, distance_comp{2,1}, .05);
if h sig = ''; else sig = 'NOT'; end
sig_str = sprintf('KS Test: Rewarded conditions are %s significantly different: h=%i p=%f', sig, h, p);
disp(sig_str);
% Also, let's throw a ranked sum test at it to test if the medians of the distributions are different
[h,p] = ranksum(distance_comp{1,1}, distance_comp{2,1}, .05);
if h sig = ''; else sig = 'NOT'; end
sig_str = sprintf('Wilcoxon: Rewarded conditions are %s significantly different: h=%i p=%f', sig, h, p);
disp(sig_str);

% Non rewarded path
yl = max([counts{1,2} counts{2,2}]);  yl = [-(yl+5) yl+5];
subplot(2,1,2); hold on;
bar(xbins, counts{1,2}, 'FaceColor', dr);
xmed = double(median(distance_comp{1,2}));
line([xmed xmed], [0 yl(1)], 'Color', dr, 'LineStyle', '--'); 
text(xmed, yl(1)+10, 1, ['median: ' num2str(xmed)], 'color', dr);
title('Trial Set 1');
set(gca, 'TickDir', 'out');
xlim(xl); ylim(yl);

bar(xbins, -counts{2,2}, 'FaceColor', dr); hold on;
xmed = double(median(distance_comp{2,2}));
line([xmed xmed], [0 -yl(1)], 'Color', dr, 'LineStyle', '--'); 
text(xmed, -yl(1)+10, ['median: ' num2str(xmed)], 'color', dr);
xlim(xl); ylim(yl);
title('Trial Set 2');
xlabel('Following Distance (px)', 'FontSize', 14);
ylabel('# Instances','FontSize', 14);
set(gca, 'TickDir', 'out');

% Plot the CDFs
plotEmpiricalCDF(distance_comp, .2, {[0 0 0],[1 0 1], [1 0 0], [1 0 0]}, plotOptions, cdf_ah);
set(get(cdf_ah,'XLabel'), 'String', 'Nose distance From Trail (px)', 'FontSize', 18);
set(get(cdf_ah,'Ylabel'), 'String', 'Cumulative Proportion of Frames', 'FontSize', 18);
set(cdf_ah, 'Tickdir', 'out', 'FontSize', 14);
axes(cdf_ah); ylim([0 1]);

% Also let's plot the n versus n+1 values against each other to see if they
% are correlated.
figure;
subplot(2,2,1);
xd = counts{1,1};
plot(xd(1:end-1), xd(2:end), 'k.'); hold on;
subplot(2,2,2);
xd = counts{2,1};
plot(xd(1:end-1), xd(2:end), 'k.'); hold on;
subplot(2,2,3);
xd = counts{1,2};
plot(xd(1:end-1), xd(2:end), 'k.'); hold on;
subplot(2,2,4);
xd = counts{2,2};
plot(xd(1:end-1), xd(2:end), 'k.'); hold on;

% Report some statistics for the distributions
for ii = 1:prod(size(counts))
    n = sum(counts{ii});
    prob = cumsum(counts{ii})./n;
    med_bin = find(prob>.5, 1, 'first');
    med = xbins(med_bin);
    low_bin = find(prob>.05, 1, 'first');
    high_bin = find(prob>.95, 1, 'first');
    disp(sprintf('For distribution %d, the median is %f, and the 5/95 percent limits are %f, %f', ...
             ii, med, xbins(low_bin), xbins(high_bin))); 
end