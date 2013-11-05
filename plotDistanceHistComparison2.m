function plotDistanceHistComparison2(dist_cell, dist_cell2, dist_thresh, plotOptions, varargin)
% function plotDistanceHists(dist_cell)
%
% Function will plot histograms of the distances from the path for 2 conditions.  They could be whatever 
% comparison.  Will plot a histogram (the pdf) on a new axis, and the cdf on either a new axis, or one
% specified by varargin. Plot options is a cell array for the linestyles (1,2) and the colors (3,4) of the cdf
% traces.
min_dist = 0; %.75;

if ~isempty(varargin)
    cdf_ah = varargin{1};
else
    figure; cdf_ah = axes; hold on;
end
if isempty(plotOptions)
   plotOptions = {'-','-','--','--'}; 
end
%defining colors
dg = [0 .8 0]; %darker green
dr = [.8 0 0]; %darker red
plotColors = {[0 0 0],dr, dg, [1 0 1]};
nfiles = [length(dist_cell); length(dist_cell2)];
rew = {dist_cell, dist_cell2};
distance_comp = cell(2,1);
for ii = 1:2
    %rewarded path
    trials = nfiles(ii);
    for jj=1:trials
        %compile distances from the trail for two conditions
        r_trials = rew{ii};
        dists = r_trials{jj};
        inci = (abs(dists) <= dist_thresh) & ~isnan(dists) & (abs(dists) >= min_dist);
        dists = dists(inci);
        distance_comp{ii} = cat(1, distance_comp{ii}, dists(:));
    end
end
    
counts = cell(2,1); %want to do a reward/distractor, early/late comparison
dx = .5;
xbins = -dist_thresh:dx:dist_thresh;
%xcenters = (dx*(0:(length(xbins)-1))) - dist_thresh + (dx/2);
counts{1} = hist(distance_comp{1}, xbins);
counts{2} = hist(distance_comp{2}, xbins); 

% plotting things
xl = [-dist_thresh dist_thresh];  %limits
yl = max([counts{1,1} counts{2,1}]); yl = [0 yl+5];

figure; ah = axes; hold on;
% Plotting PDF of first data set
%bar(xbins, counts{1}, 'FaceColor', dg); hold on;
xmed = double(median(distance_comp{1}));
line([xmed xmed], [0 yl(2)], 'Color',plotColors{1}, 'LineStyle', '--'); 
text(xmed, yl(2), ['median: ' num2str(xmed)], 'Color', plotColors{1});
line(xbins, counts{1}, 'Color', plotColors{1}, 'LineStyle', '-', 'LineWidth', 2);

%bar(xbins, counts{2}, 'FaceColor', dg); hold on;
xmed = double(median(distance_comp{2}));
line([xmed xmed], [0 yl(2)], 'Color', plotColors{2}, 'LineStyle', '--'); 
text(xmed, yl(2)/2, ['median: ' num2str(xmed)], 'color', plotColors{2});
line(xbins, counts{2}, 'Color', plotColors{2}, 'LineStyle', '-', 'LineWidth', 2);

% KS test the distributions
if ~isempty(distance_comp{1}) && ~isempty(distance_comp{2})
    [h,p,ks2stat] = kstest2(distance_comp{1}, distance_comp{2}, .05);
    if h sig = ''; else sig = 'NOT'; end
    sig_str = sprintf('KS Test: Conditions are %s significantly different: h=%i p=%f', sig, h, p);
    disp(sig_str);
    % Also, let's throw a ranked sum test at it to test if the medians of the distributions are different
    [h,p] = ranksum(distance_comp{1}, distance_comp{2}, .05);
    if h sig = ''; else sig = 'NOT'; end
    sig_str = sprintf('Wilcoxon: Conditions are %s significantly different: h=%i p=%f', sig, h, p);
    disp(sig_str);
end

for ii=1:length(distance_comp)
    if ~isempty(distance_comp{ii})
        plotEmpiricalCDF(distance_comp(ii), .2, plotColors(ii), plotOptions(ii), cdf_ah);
    end
end
set(get(cdf_ah,'XLabel'), 'String', 'Nose distance From Trail (px)', 'FontSize', 18);
set(get(cdf_ah,'Ylabel'), 'String', 'Cumulative Proporation of Frames', 'FontSize', 18);
set(cdf_ah, 'Tickdir', 'out', 'FontSize', 14);
axes(cdf_ah); ylim([0 1]);

% Non rewarded path
% yl = max([counts{1,2} counts{2,2}]);  yl = [-(yl+5) yl+5];
% subplot(2,1,2); hold on;
% bar(xbins, counts{1,2}, 'FaceColor', dr);
% xmed = double(median(distance_comp{1,2}));
% line([xmed xmed], [0 yl(1)], 'Color', dr, 'LineStyle', '--'); 
% text(xmed, yl(1)+10, 1, ['median: ' num2str(xmed)], 'color', dr);
% title('Trial Set 1');
% set(gca, 'TickDir', 'out');
% xlim(xl); ylim(yl);
% 
% bar(xbins, -counts{2,2}, 'FaceColor', dr); hold on;
% xmed = double(median(distance_comp{2,2}));
% line([xmed xmed], [0 -yl(1)], 'Color', dr, 'LineStyle', '--'); 
% text(xmed, -yl(1)+10, ['median: ' num2str(xmed)], 'color', dr);
% xlim(xl); ylim(yl);
% title('Trial Set 2');
% xlabel('Following Distance (px)', 'FontSize', 14);
% ylabel('# Instances','FontSize', 14);
% set(gca, 'TickDir', 'out');

% % Plot the CDFs
% plotEmpiricalCDF(distance_comp, .2, {[0 0 0],[1 0 1], [1 0 0], [1 0 0]}, plotOptions, cdf_ah);
% set(get(cdf_ah,'XLabel'), 'String', 'Nose distance From Trail (px)', 'FontSize', 18);
% set(get(cdf_ah,'Ylabel'), 'String', 'Cumulative Proporation of Frames', 'FontSize', 18);
% set(cdf_ah, 'Tickdir', 'out', 'FontSize', 14);
% axes(cdf_ah); ylim([0 1]);

% Also let's plot the n versus n+1 values against each other to see if they
% are correlated.
% figure;
% subplot(2,2,1);
% xd = counts{1,1};
% plot(xd(1:end-1), xd(2:end), 'k.'); hold on;
% subplot(2,2,2);
% xd = counts{2,1};
% plot(xd(1:end-1), xd(2:end), 'k.'); hold on;
% subplot(2,2,3);
% xd = counts{1,2};
% plot(xd(1:end-1), xd(2:end), 'k.'); hold on;
% subplot(2,2,4);
% xd = counts{2,2};
% plot(xd(1:end-1), xd(2:end), 'k.'); hold on;
