% fallCumDistPlots     
% --------------------
% The point of this is to analyze the mouse position relative to the trail
% as a time-series, using models that can then be used to draw statistical
% conclusions from it.  The bias towards one side of the trail or the
% other, is a goal.  To perform the analysis, we loop through the
% individual mice and each trial to pick out when they were following,
% concatenate each segment of following together, and use that as our
% time-series data.

% Parameters
nMice = length(perMouseData);
nRows = ceil(sqrt(nMice));
distThresh = 20;
mm_conv = .862; %mm/px conversion
distThresh = distThresh*mm_conv;
segLenThresh = 40; % n samples in a row minimum for this analysis
AR_width = 4;
armean = @(fitp) fitp(1) ./ (1-sum(fitp(2:end))); %mean of the AR model
anonARmean = @(data) armean(fitARmodel(data, AR_width));
nboot = 200; %number of bootstrap resamples
segL = 100; %resampling segment length for bootstrapping
segSep = 100;


% Initialize some vars
fitp = cell(nMice, 4); fitp_ci = cell(nMice, 4); %parameters of model fit
armean = cell(nMice, 4); armean_ci = cell(nMice, 4); %parameters of model fit
ar_means = NaN*zeros(nMice,4);
e = cell(nMice, 4); fit = cell(nMice, 4);
med_vel = NaN*zeros(nMice, 4); %median velocities - want to check these
acorrFigh = figure;
segLenFigh = figure;
residFigh = figure;
cumfigh = figure;

for jj = 1:nMice
    figure(cumfigh);
    subplot(nRows, nRows, jj); hold on;
    ts = []; trialCount = []; segLen = [];
    mouseDists = perMouseData(jj).rew_dists_from_trail;
    for ii = 1:length(mouseDists)
        % Plot some markers for epoch beginnings
        xp = length(ts); 
        if ~isempty(ctl_trials{jj}) && ctl_trials{jj}(1) == ii
            plot([xp xp], [-10 10].*1e4, 'k--');
        end
        if ~isempty(ctl2_trials{jj})
            if ctl2_trials{jj}(1) == ii
                plot([xp xp], [-10 10].*1e4, 'k--');
            end
        end
        if ~isempty(occr_trials{jj})
            if occr_trials{jj}(1) ==  ii
                plot([xp xp], [-10 10].*1e4, 'b--');
            end
        end
        if ~isempty(occl_trials{jj})
            if occl_trials{jj}(1) == ii
                plot([xp xp], [-10 10].*1e4, 'r--');
            end
        end
        % want to select for following segments, but we also want to make
        % sure we exclude the really short segments because concatenating too
        % many short segments could hurt our efforts to model the behavior.  
        seli = find(abs(mouseDists{ii}) <= distThresh);
        segs = findContinuousSegments(seli); % n x 2 matrix with start and end of each seg
        segLenTmp = diff(segs, [], 2) + 1; 
        long = segLenTmp >= segLenThresh;
        segs = segs(long,:);
        segi = [];
        for kk = 1:size(segs,1)
            segi = cat(1, segi, (segs(kk,1):segs(kk,2))' ); 
        end
        segLen = cat(1, segLen, segLenTmp(long));
        following = mouseDists{ii}(segi);
        %following = mouseDists{ii}(abs(mouseDists{ii}) <= distThresh); % easier to just do it again, logically
        ts = cat(1, ts, following);
        trialCount = cat(1, trialCount, ii*ones(length(following),1));
    end
    % Let's do some Timeseries object stuff with this - the advantage of
    % this is that it linearly interprets missing data (NaN).
    posTS = timeseries(ts, 1:length(ts), 'Name', 'Distance From Trail'); 
    interpTS = posTS.Data;
    ts = interpTS;
    %ts = ts(~isnan(ts));
    cum_ts = cumsum(ts);
    %plot the data, fix up the axes
    min_ts = min(cum_ts); max_ts = max(cum_ts); pad = 1000;
    ylim([min_ts - pad,  max_ts+pad]);
    plot(cum_ts, 'LineWidth', 1, 'Color', 'k');
    set(gca, 'Ydir', 'reverse');
    ylabel('Cumulative Distance From Trail (mm)');
    xlabel('Behavioral Time (frames)');
    
    figure(segLenFigh);
    subplot(nRows, nRows, jj); hold on;
    hist(segLen, 100);
%     posTS = timeseries(ts, 1:length(ts), 'Name', 'Distance From Trail'); 
%     posTS.DataInfo.Units = 'mm';
%     tscol = tscollection(posTS, 'Name', mouse_names{jj});
%     cumTS = timeseries(cumsum(posTS.Data), 1:length(ts), 'Name', 'Cumulative Distance');
%     tscol = addts(tscol, cumTS);
%     if ~isempty(ctl_trials{jj})
%         epochi = find(trialCount >= ctl_trials{jj}(1) & trialCount <= ctl_trials{jj}(end));
%         data = posTS.Data(epochi);
%         sum(isnan(data));
%         [fitp, fitp_ci, e, fit] = fitARmodel(data, 5);
%     end
    
    % Want to do some linear fits to the various training epochs
    if ~isempty(ctl_trials{jj})
        fprintf('Control: ');
        epochi = find(trialCount >= ctl_trials{jj}(1) & trialCount <= ctl_trials{jj}(end));
        %linear trend fit of cumulative
        fitx = epochi - epochi(1);
        [p, fity, se, ci] = linear_fit_stats(fitx, cum_ts(epochi)); %my linear fitting with CIs added (note, assumes IID)
        %plot(epochi, fity,'Color',[.5 .5 .5], 'LineWidth',1);
        slopes(jj,1) = p(1); slope_ci(jj,1) = ci(1);
        
        % fitting autoregressive (AR with mean) model to timeseries
        [fitp{jj,1}, fitp_ci{jj,1}, e{jj,1}, fit{jj,1}] = fitARmodel(ts(epochi), AR_width);
        param_sum = sum(fitp{jj,1}(2:end));  fprintf('Unit root of fit is %f\n', param_sum);
        ar_means(jj,1) =  fitp{jj,1}(1) ./ (1-sum(fitp{jj,1}(2:end))); %mean of the AR model
        % putting confidence intervals on the mean using moving block bootstrapping
        [armean{jj,1}, armean_ci{jj,1}] = mbBootstrap(nboot, anonARmean, ts(epochi), segL, segSep);
        epochy = cum_ts(epochi);
        predy =  ((1:length(epochy)) * ar_means(jj,1))'; 
        p = polyfit(fitx, epochy-predy, 0);
        figure(cumfigh);
        plot(epochi, predy+p,'Color',[.5 .5 .5], 'LineWidth',1);
        figure(residFigh);
        subplot(nMice, 4, 4*(jj-1)+1); hold on;
        plot(e{jj,1}(1:end-1), e{jj,1}(2:end),'k.');
        % Plot the autocorrelation function of this dataset
        figure(acorrFigh);
        subplot(nRows, nRows, jj); hold on;
        [acorr, lags] = xcorr(ts(epochi), 30, 'unbiased');
        plot(lags, acorr, 'k');
        med_vel(jj,1) = median(abs(diff(ts(epochi))));
        %plot(epochi((AR_width+2):end), cum_ts(epochi(1)+AR_width-1)+cumsum(fit{jj,1}), '--', 'Color',[.5 .5 .5]);
    end
    if ~isempty(ctl2_trials{jj})
        fprintf('Control2: ');
        epochi = find(trialCount >= ctl2_trials{jj}(1) & trialCount <= ctl2_trials{jj}(end));
        fitx = epochi - epochi(1);
        [p, fity, se, ci] = linear_fit_stats(fitx, cum_ts(epochi));
        %plot(epochi, fity,'Color',[.5 .5 .5], 'LineWidth',1);
        slopes(jj,2) = p(1); slope_ci(jj,2) = ci(1);
        
        % fitting autoregressive (AR with mean) model to timeseries
        [fitp{jj,2}, fitp_ci{jj,2}, e{jj,2}, fit{jj,2}] = fitARmodel(ts(epochi), AR_width);
        param_sum = sum(fitp{jj,2}(2:end));  fprintf('Unit root of fit is %f\n', param_sum);
        ar_means(jj,2) =  fitp{jj,2}(1) ./ (1-sum(fitp{jj,2}(2:end))); %mean of the AR model
        [armean{jj,2}, armean_ci{jj,2}] = mbBootstrap(nboot, anonARmean, ts(epochi), segL, segSep);
        epochy = cum_ts(epochi);
        predy =  ((1:length(epochy)) * ar_means(jj,2))'; 
        p = polyfit(fitx, epochy-predy, 0);
        figure(cumfigh);
        plot(epochi, predy+p,'Color',[.5 .5 .5], 'LineWidth',1);
        %plot(epochi((AR_width+2):end), cum_ts(epochi(1)+AR_width-1)+cumsum(fit{jj,2}), '--', 'Color',[.5 .5 .5]);
        figure(residFigh);
        subplot(nMice, 4, 4*(jj-1)+2); hold on;
        plot(e{jj,2}(1:end-1), e{jj,2}(2:end),'.', 'Color', [.5 .5 .5]);
        % Plot the autocorrelation function of this dataset
        figure(acorrFigh); hold on;
        [acorr, lags] = xcorr(ts(epochi), 30, 'unbiased');
        plot(lags, acorr, 'Color', [.5 .5 .5]);
        med_vel(jj,2) = median(abs(diff(ts(epochi))));
    else
        slopes(jj,2) = NaN; slope_ci(jj,2) = NaN;
    end
    if ~isempty(occr_trials{jj})
        fprintf('Right Occlusion: ');
        epochi = find(trialCount >= occr_trials{jj}(1) & trialCount <= occr_trials{jj}(end));
        fitx = epochi - epochi(1);
        [p, fity, se, ci] = linear_fit_stats(fitx, cum_ts(epochi));
        %plot(epochi, fity,'b', 'LineWidth',1);
        slopes(jj,3) = p(1); slope_ci(jj,3) = ci(1);
        
        % fitting autoregressive (AR with mean) model to timeseries
        [fitp{jj,3}, fitp_ci{jj,3}, e{jj,3}, fit{jj,3}] = fitARmodel(ts(epochi), AR_width);
        param_sum = sum(fitp{jj,3}(2:end));  fprintf('Unit root of fit is %f\n', param_sum);
        ar_means(jj,3) =  fitp{jj,3}(1) ./ (1-sum(fitp{jj,3}(2:end))); %mean of the AR model
        [armean{jj,3}, armean_ci{jj,3}] = mbBootstrap(nboot, anonARmean, ts(epochi), segL, segSep);
        epochy = cum_ts(epochi);
        predy =  ((1:length(epochy)) * ar_means(jj,3))'; 
        p = polyfit(fitx, epochy-predy, 0);
        figure(cumfigh);
        plot(epochi, predy+p,'Color','b', 'LineWidth',1);
        %plot(epochi((AR_width+2):end), cum_ts(epochi(1)+AR_width-1)+cumsum(fit{jj,3}), '--', 'Color',[.5 .5 .5]);
        figure(residFigh);
        subplot(nMice, 4, 4*(jj-1)+3); hold on;
        plot(e{jj,3}(1:end-1), e{jj,3}(2:end),'b.');
        % Plot the autocorrelation function of this dataset
        figure(acorrFigh); hold on;
        [acorr, lags] = xcorr(ts(epochi), 30, 'unbiased');
        plot(lags, acorr, 'b');
        med_vel(jj,3) = median(abs(diff(ts(epochi))));
    else
        slopes(jj,3) = NaN; slope_ci(jj,3) = NaN;
    end
    if ~isempty(occl_trials{jj})
        fprintf('Left Occlusion: ');
        epochi = find(trialCount >= occl_trials{jj}(1) & trialCount <= occl_trials{jj}(end));
        fitx = epochi - epochi(1);
        [p, fity, se, ci] = linear_fit_stats(fitx, cum_ts(epochi));
        %plot(epochi, fity,'r', 'LineWidth',1);
        slopes(jj,4) = p(1); slope_ci(jj,4) = ci(1);
        
         % fitting autoregressive (AR with mean) model to timeseries
        [fitp{jj,4}, fitp_ci{jj,4}, e{jj,4}, fit{jj,4}] = fitARmodel(ts(epochi), AR_width);
        param_sum = sum(fitp{jj,4}(2:end));  fprintf('Unit root of fit is %f\n', param_sum);
        ar_means(jj,4) =  fitp{jj,4}(1) ./ (1-sum(fitp{jj,4}(2:end))); %mean of the AR model
        [armean{jj,4}, armean_ci{jj,4}] = mbBootstrap(nboot, anonARmean, ts(epochi), segL, segSep);
        epochy = cum_ts(epochi);
        predy =  ((1:length(epochy)) * ar_means(jj,4))'; 
        p = polyfit(fitx, epochy-predy, 0);
        figure(cumfigh);
        plot(epochi, predy+p,'Color','r', 'LineWidth',1);
        %plot(epochi((AR_width+2):end), cum_ts(epochi(1)+AR_width-1)+cumsum(fit{jj,4}), '--', 'Color',[.5 .5 .5]);
        figure(residFigh);
        subplot(nMice, 4, 4*(jj-1)+4); hold on;
        plot(e{jj,4}(1:end-1), e{jj,4}(2:end),'r.');
        % Plot the autocorrelation function of this dataset
        figure(acorrFigh); hold on;
        [acorr, lags] = xcorr(ts(epochi), 30, 'unbiased');
        plot(lags, acorr, 'r');
        med_vel(jj,4) = median(abs(diff(ts(epochi))));
    else
        slopes(jj,4) = NaN; slope_ci(jj,4) = NaN;
    end
    % Check if the confidence intervals overlap
    
end
slopes
slope_ci

ar_slopes = NaN*zeros(size(fitp));
slope_lb =  NaN*zeros(size(fitp));
slope_ub =  NaN*zeros(size(fitp));
ar_means = NaN*zeros(size(fitp));
for ii = 1:numel(fitp)
    temp = fitp{ii};
    tci = fitp_ci{ii};
    if ~isempty(temp)
        ar_slopes(ii) = temp(1);
        slope_lb(ii) = tci(1,1);
        slope_ub(ii) = tci(1,2);
        ar_means(ii) = ar_slopes(ii) ./ (1-sum(temp(2:end)));
    end
end
