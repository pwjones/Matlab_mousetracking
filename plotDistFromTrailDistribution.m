function fitp = plotDistFromTrailDistribution(dist_cell, dist_thresh, plotOptions, varargin)
% function plotDistanceHists(dist_cell)
%
% Function will plot histograms of the distances from the path.
% Will plot a histogram (the pdf) on a new axis, and the cdf on either a new axis, or one
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
nfiles = length(dist_cell);
rew = {dist_cell};
all_dists = [];

ii =1;
trials = nfiles(ii);
for jj=1:trials
    %compile distances from the trail for two conditions
    r_trials = rew{ii};
    dists = r_trials{jj};
    inci = (abs(dists) <= dist_thresh) & ~isnan(dists) & (abs(dists) >= min_dist);
    dists = dists(inci);
    all_dists = cat(1, all_dists, dists(:));
end
    
nbins = 10*dist_thresh;
mlePDF = @cauchyPDF2;
mleCDF = @cauchyCDF2;

xbins = linspace(0,dist_thresh,floor(nbins/2)+1); %might be an easier way - symmetric vector around 0
rev = -xbins(2:end); xbins = [rev(end:-1:1), xbins];
x0 = [nanmean(all_dists), 1]; 
% make a PDF of the distribution of distances
counts = hist(all_dists, xbins);
count_dist = counts./nansum(counts);

% use mle rather than lsqcurvefit for fitting a distribution - it's supposedly the right 
mle_fit_params = mle(all_dists, 'pdf', mlePDF, 'start', double(x0), 'cdf', mleCDF);
[f, xcdf] = plotEmpiricalCDF({all_dists}, .2, plotColors(ii), plotOptions(ii), cdf_ah);
gfit = mleCDF(xbins, mle_fit_params(1), mle_fit_params(2));
line('parent', cdf_ah,'xdata',xbins, 'ydata', gfit, 'Linestyle', '--', 'color','r');
disp(sprintf('Full Mean: %.3f   SD: %.3f', nanmean(all_dists), nanstd(all_dists)));
tail_dists = all_dists(all_dists ~= 0);
disp(sprintf('Tail Data Mean: %.3f   SD: %.3f', nanmean(tail_dists), nanstd(tail_dists)));
disp(['\nMLE Fitted Cauchy dist has parameters: ' num2str(mle_fit_params(1)) '  And width: ' num2str(mle_fit_params(2))]);


% Plotting limits
xl = [-dist_thresh dist_thresh];  
yl = counts./sum(counts); 
yl = [0 yl+(yl/5)];
% Plotting PDFs 
figure; ah = axes; hold on;
for ii =1:1
    % The PDF 
    line(xbins, count_dist, 'Color', plotColors{ii}, 'LineStyle', '-', 'LineWidth', 2);
    xmed = double(median(all_dists)); 
    line([xmed xmed], [0 yl(2)], 'Color',plotColors{ii}, 'LineStyle', '--'); 
    text(xmed, yl(2), ['median: ' num2str(xmed)], 'Color', plotColors{ii});
    % The fitted function
    gfit = mlePDF(xbins, mle_fit_params(1), mle_fit_params(2));
    line(xbins, gfit./sum(gfit), 'Linestyle', '--', 'color','r');
end

set(get(cdf_ah,'XLabel'), 'String', 'Nose distance From Trail (px)', 'FontSize', 18);
set(get(cdf_ah,'Ylabel'), 'String', 'Cumulative Proporation of Frames', 'FontSize', 18);
set(cdf_ah, 'Tickdir', 'out', 'FontSize', 14);
axes(cdf_ah); ylim([0 1]); xlim(xl);

fitp = mle_fit_params;
