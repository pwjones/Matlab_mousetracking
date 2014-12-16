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
distance_comp = cell(1,1);

ii =1;
trials = nfiles(ii);
for jj=1:trials
    %compile distances from the trail for two conditions
    r_trials = rew{ii};
    dists = r_trials{jj};
    inci = (abs(dists) <= dist_thresh) & ~isnan(dists) & (abs(dists) >= min_dist);
    dists = dists(inci);
    distance_comp{ii} = cat(1, distance_comp{ii}, dists(:));
end
    
counts = cell(1,1); %want to do a reward/distractor, early/late comparison
nbins = 10*dist_thresh;
dx = .5;
%anonPDF = @(x,xdata)cauchyPDF2(xdata, x(1), x(2));
%anonPDF = @(x,xdata)nctpdf(xdata, x(2), x(1)); % the mean and width parameters are switched
mlePDF = @cauchyPDF2;
mleCDF = @cauchyCDF2;
%mlePDF = @normpdf;
%mleCDF = @normcdf;
%xbins = linspace(-dist_thresh, dist_thresh, nbins);
xbins = linspace(0,dist_thresh,floor(nbins/2)+1); %might be an easier way - symmetric vector around 0
rev = -xbins(2:end); xbins = [rev(end:-1:1), xbins];
x0 = [nanmean(distance_comp{1}), 1]; 
opts =  optimset('MaxFunEvals', 5000, 'TolFun', 1e-15, 'MaxIter', 10000, 'Diagnostics', 'on');
%xcenters = (dx*(0:(length(xbins)-1))) - dist_thresh + (dx/2);
counts = hist(distance_comp{1}, xbins);
%zeroi = find(xbins > 0, 1,'first'); %jsut a test for what is messing up fits
%max_count = max(counts);
%counts((zeroi-1):zeroi) = max_count;
count_dist = counts./nansum(counts);
fit_params = [0 1]; fit_paramsCDF = [0,1];
%[fit_params, resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(anonPDF, double(x0), xbins, count_dist,[],[10 10],opts);
mle_fit_params = mle(distance_comp{1}, 'pdf', mlePDF, 'start', double(x0), 'cdf', mleCDF);
%[mle_fit_params(1), mle_fit_params(2)] = normfit(distance_comp{1});
%anonCDF = @(x,xdata)cauchyCDF2(xdata, x(1), x(2));
%anonCDF = @(x,xdata)nctcdf(xdata, x(2), x(1));
[f, xcdf] = plotEmpiricalCDF(distance_comp(ii), .2, plotColors(ii), plotOptions(ii), cdf_ah);
%[fit_paramsCDF, resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(anonCDF, double(x0), double(xcdf), double(f),[],[10 30],opts);
% use mle rather than lsqcurvefit
%gfit = anonCDF(fit_params, xbins);
%line('parent', cdf_ah,'xdata',xbins, 'ydata', gfit, 'Linestyle', '--', 'color','r');
gfit = mleCDF(xbins, mle_fit_params(1), mle_fit_params(2));
line('parent', cdf_ah,'xdata',xbins, 'ydata', gfit, 'Linestyle', '--', 'color','r');
   
%disp(['Fitted Anon dist PDF has parameters: ' num2str(fit_params(1)) '  And width: ' num2str(fit_params(2))]);
%disp(['Fitted Anon dist CDF has parameters: ' num2str(fit_paramsCDF(1)) '  And width: ' num2str(fit_paramsCDF(2))]);
disp(['Fitted MLE dist CDF has parameters: ' num2str(mle_fit_params(1)) '  And width: ' num2str(mle_fit_params(2))]);

% Plotting limits
xl = [-dist_thresh dist_thresh];  
yl = counts./sum(counts); 
yl = [0 yl+(yl/5)];
% Plotting PDFs 
figure; ah = axes; hold on;
for ii =1:1
    % The PDF 
    line(xbins, count_dist, 'Color', plotColors{ii}, 'LineStyle', '-', 'LineWidth', 2);
    xmed = double(median(distance_comp{ii})); 
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
