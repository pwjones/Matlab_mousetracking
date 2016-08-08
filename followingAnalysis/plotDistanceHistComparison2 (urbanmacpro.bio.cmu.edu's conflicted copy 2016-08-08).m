function plotDistanceHistComparison2(dist_cell, dist_cell2, dist_thresh, plotOptions, varargin)
% function plotDistanceHists(dist_cell)
%
% Function will plot histograms of the distances from the path for 2 conditions.  They could be whatever 
% comparison.  Will plot a histogram (the pdf) on a new axis, and the cdf on either a new axis, or one
% specified by varargin. Plot options is a cell array for the linestyles (1,2) and the colors (3,4) of the cdf
% traces.
min_dist = 0; %.75;
%defining colors
dg = [0 .8 0]; %darker green
dr = [.8 0 0]; %darker red
plotColors = {[0 0 0],dr, dg, [1 0 1]};

%Checking variable arguments
if ~isempty(varargin)
    cdf_ah = varargin{1};
else
    figure; cdf_ah = axes; hold on;
end
if nargin > 5 %2 or more variable args
    plotColors = varargin{2};
end
if isempty(plotOptions)
   plotOptions = {'-','-','--','--'}; 
end

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
nbins = 60;
dx = .5;
xbins = linspace(-dist_thresh, dist_thresh, nbins);
%xcenters = (dx*(0:(length(xbins)-1))) - dist_thresh + (dx/2);
for ii=1:2
    counts{ii} = hist(distance_comp{ii}, xbins);
    [muhat(ii), sigmahat(ii)] = normfit(distance_comp{ii});
    mu(ii) = nanmean(distance_comp{ii});
    sigma(ii) = nanstd(distance_comp{ii});
    med(ii) = nanmedian(distance_comp{ii});
end
fprintf('Data distributions have Medians: %s  Means: %s, STDs: %s\n', num2str(med), num2str(mu), num2str(sigma));
disp(['Fitted Gaussians have means: ' num2str(muhat) '  And stds: ' num2str(sigmahat)]);

% plotting things
xl = [-dist_thresh dist_thresh];  %limits
yl = max([counts{1,1}./sum(counts{1,1}), counts{2,1}./sum(counts{2,1})]); yl = [0 yl+(yl/5)];

figure; ah = axes; hold on;
% Plotting PDFs 
for ii =1:2
    xmed = double(median(distance_comp{ii}));
    xmean = mu(ii);
    line([xmean xmean], [0 yl(2)], 'Color',plotColors{ii}, 'LineStyle', '--'); 
    text(xmean, yl(2), ['Mean: ' num2str(xmean)], 'Color', plotColors{ii});
    line([xmed xmed], [0 yl(2)], 'Color',plotColors{ii}, 'LineStyle', '-'); 
    text(xmed, yl(2), ['Median: ' num2str(xmed)], 'Color', plotColors{ii});
    line(xbins, counts{ii}./sum(counts{ii}), 'Color', plotColors{ii}, 'LineStyle', '-', 'LineWidth', 2);
    gfit = normpdf(xbins, muhat(ii), sigmahat(ii));
    line(xbins, gfit./sum(gfit), 'Linestyle', ':', 'color',plotColors{ii});
end

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
        gfit = normcdf(xbins, muhat(ii), sigmahat(ii));
        line('parent', cdf_ah,'xdata',xbins, 'ydata', gfit./gfit(end), 'Linestyle', '--', 'color',plotColors{ii});
    end
end
set(get(cdf_ah,'XLabel'), 'String', 'Nose distance From Trail (px)', 'FontSize', 18);
set(get(cdf_ah,'Ylabel'), 'String', 'Cumulative Proporation of Frames', 'FontSize', 18);
set(cdf_ah, 'Tickdir', 'out', 'FontSize', 14);
axes(cdf_ah); ylim([0 1]); xlim(xl);
