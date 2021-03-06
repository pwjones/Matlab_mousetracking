% fallCumDistPlots     
% --------------------
% The point of this is to analyze the mouse position relative to the trail
% as a time-series, using models that can then be used to draw statistical
% conclusions from it.  The bias towards one side of the trail or the
% other, is a goal.  To perform the analysis, we loop through the
% individual mice and each trial to pick out when they were following,
% concatenate each segment of following together, and use that as our
% time-series data.

% What to plot 
plotAcorr = 0; 
plotSegHist = 0;
plotResid = 0;
% Parameters
nMice = length(perMouseData);
nRows = ceil(sqrt(nMice));
distThresh = 20;
mm_conv = .862; %mm/px conversion
distThresh = distThresh*mm_conv;
segLenThresh = 40; % n samples in a row minimum for this analysis
corrWidth = 40;
AR_width = 3;
MA_width = 4;
nboot = 500; %number of bootstrap resamples
segL = 40; %resampling segment length for bootstrapping
segSep = 50;
% anonymous functions used to bootstrap the means of the fitted AR models
armean = @(fitp) fitp(1) ./ (1-sum(fitp(2:end))); %mean of the AR model from fitted parameters
anonARmean = @(data) armean(fitARmodel(data, AR_width)); % returns the mean rather than all params, can call in loop easily


% Initialize some vars
fitp = cell(nMice, 4); fitp_ci = cell(nMice, 4); %parameters of model fit
armean = NaN*zeros(nMice, 4); armean_ci = NaN*zeros(nMice, 4, 2); %parameters of model fit
allmeans = NaN*zeros(nMice, 4, nboot);
dwstat = NaN*zeros(nMice,4);
ar_means = NaN*zeros(nMice,4);
e = cell(nMice, 4); fit = cell(nMice, 4);
med_vel = NaN*zeros(nMice, 4); %median velocities - want to check these
if (plotAcorr) acorrFigh = figure; end
if (plotSegHist) segLenFigh = figure; end
if (plotResid) residFigh = figure; end
spectfigh = figure; axes;
cumfigh = figure;



for jj = 1:nMice
    figure(cumfigh);
    subplot(nRows, nRows, jj); hold on;
    ts = []; trialCount = []; segLen = [];
    followingSegs = {};
    mouseDists = perMouseData(jj).rew_dists_from_trail; % this is the distance from each individual trial/video
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
        % want to select for following segments, but we also want to make sure we exclude the really 
        % short segments because we don't want to pollute the timeseries with too many introduced discontinuous pieces.  
        seli = find(abs(mouseDists{ii}) <= distThresh);
        segs = findContinuousSegments(seli); % n x 2 matrix with start and end of each seg
        segLenTmp = diff(segs, [], 2) + 1; 
        long = segLenTmp >= segLenThresh;
        segs = segs(long,:);
        segi = [];
        for kk = 1:size(segs,1) % make the selection vector
            segi = cat(1, segi, (segs(kk,1):segs(kk,2))' ); 
        end
        for kk = 1:size(segs, 1) % we also want to make a cell array of the inidividual segments
           si =  segs(kk,1):segs(kk,2);
           followingSegs{end+1} = mouseDists{ii}(si);
        end
        segLen = cat(1, segLen, segLenTmp(long));
        following = mouseDists{ii}(segi);
        ts = cat(1, ts, following);
        trialCount = cat(1, trialCount, ii*ones(length(following),1)); %vector giving the trail that each segment comes from
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
    set(gca, 'Ydir', 'reverse', 'TickDir', 'out');
    ylabel('Cumulative Distance From Trail (mm)');
    xlabel('Behavioral Time (frames)');
    
    % Spectral analysis of the following nose positions in order to look at
    % scanning of the trail
    distMat = makePaddedMatFromCell(followingSegs);
    distMat_filt = gaussianFilter(distMat, 3);
    [f,S] = noseSpectAnal2(distMat_filt,1/40);
    figure(spectfigh); hold on; plot(f,S);
    
    if (plotSegHist)
        figure(segLenFigh);
        subplot(nRows, nRows, jj); hold on;
        hist(segLen, 100);
    end
    
    % Want to do some linear fits to the various training epochs
    if ~isempty(ctl_trials{jj})
        fprintf('Control: ');
        % select indices based on the trails numbers that they are pulled from to get the right experimental condition 
        epochi = find(trialCount >= ctl_trials{jj}(1) & trialCount <= ctl_trials{jj}(end));
        
        
%         [f,S, tah, fah] = noseSpectAnal2(ts(epochi),1/40);
%         figure(spectfigh); hold on; plot(f,S);
%         tsbreaks = cumsum(segLen);
%         if isempty(tah)
%             figure; tah = axes; hold on;
%         end
%         axes(tah); hold on; 
%         plot(tah, tsbreaks./40, zeros(size(tsbreaks)), 'r.');
        
        %linear trend fit of cumulative
        fitx = epochi - epochi(1);
        [p, fity, se, ci] = linear_fit_stats(fitx, cum_ts(epochi)); %my linear fitting with CIs added (note, assumes IID)
        %plot(epochi, fity,'Color',[.5 .5 .5], 'LineWidth',1);
        slopes(jj,1) = p(1); slope_ci(jj,1) = ci(1);
        ctrl_ts = ts(epochi);
        y = ts(epochi);
        
        % fitting autoregressive (AR with mean) model to timeseries
        [fitp{jj,1}, fitp_ci{jj,1}, e{jj,1}, fit{jj,1}, dwstat(jj,1)] = fitARmodel(ts(epochi), AR_width);
        param_sum = sum(fitp{jj,1}(2:end));  fprintf('Unit root of fit is %f\n', param_sum);
        ar_means(jj,1) =  fitp{jj,1}(1) ./ (1-sum(fitp{jj,1}(2:end))); %mean of the AR model
        
        % putting confidence intervals on the mean using moving block bootstrapping
        [armean(jj,1), armean_ci(jj,1,:), allmeans(jj,1,:)] = mbBootstrap(nboot, anonARmean, ts(epochi), segL, segSep);
        %[armean(jj,1), armean_ci(jj,1,:)] = mbBootstrap(nboot, @fitCumTsSlope, ts(epochi), segL, segSep);
        % Using stationary bootstrap - just going to get means directly from bootstrapped data
        Bstar = opt_block_length_REV_dec07(ts(epochi)); % gives the optimal length for bootstrap blocks
        [ystar, mustar] = stationaryBootstrap(nboot, ceil(Bstar(1)), y);
        %[armean(jj,1), armean_ci(jj,1,:)] = eCI(mustar(:), .95);
        %allmeans(jj,1,:) = mustar;
        
        epochy = cum_ts(epochi);
        predy =  ((1:length(epochy)) * ar_means(jj,1))'; 
        p = polyfit(fitx, epochy-predy, 0);
        figure(cumfigh);
        plot(epochi, predy+p,'Color',[.5 .5 .5], 'LineWidth',1);
        figure(residFigh);
        subplot(nMice, 4, 4*(jj-1)+1); hold on;
        plot(e{jj,1}(1:end-1), e{jj,1}(2:end),'k.');
        % Plot the autocorrelation function of this dataset
        if plotAcorr
            figure(acorrFigh);
            subplot(nRows, nRows, jj); hold on;
            [acorr, lags] = xcorr(ts(epochi), corrWidth, 'unbiased');
            plot(lags, acorr, 'k');
        end
        med_vel(jj,1) = median(abs(diff(ts(epochi))));
        %plot(epochi((AR_width+2):end), cum_ts(epochi(1)+AR_width-1)+cumsum(fit{jj,1}), '--', 'Color',[.5 .5 .5]);
    end
    if ~isempty(ctl2_trials{jj})
        fprintf('Control2: ');
        epochi = find(trialCount >= ctl2_trials{jj}(1) & trialCount <= ctl2_trials{jj}(end));
        y = ts(epochi);
        %y = cat(1, ctrl_ts(:), y(:)); % add the two together.
        fitx = epochi - epochi(1);
        [p, fity, se, ci] = linear_fit_stats(fitx, cum_ts(epochi));
        %plot(epochi, fity,'Color',[.5 .5 .5], 'LineWidth',1);
        slopes(jj,2) = p(1); slope_ci(jj,2) = ci(1);
        
        % fitting autoregressive (AR with mean) model to timeseries
        [fitp{jj,2}, fitp_ci{jj,2}, e{jj,2}, fit{jj,2}, dwstat(jj,2)] = fitARmodel(ts(epochi), AR_width);
        param_sum = sum(fitp{jj,2}(2:end));  fprintf('Unit root of fit is %f\n', param_sum);
        ar_means(jj,2) =  fitp{jj,2}(1) ./ (1-sum(fitp{jj,2}(2:end))); %mean of the AR model
        [armean(jj,2), armean_ci(jj,2,:), allmeans(jj,2,:)] = mbBootstrap(nboot, anonARmean, ts(epochi), segL, segSep);
        %[armean(jj,2), armean_ci(jj,2,:)] = mbBootstrap(nboot, @fitCumTsSlope, ts(epochi), segL, segSep);
        % Using stationary bootstrap - just going to get means directly from bootstrapped data
        Bstar = opt_block_length_REV_dec07(ts(epochi)); % gives the optimal length for bootstrap blocks
        [ystar, mustar] = stationaryBootstrap(nboot, ceil(Bstar(1)), y);
        %[armean(jj,2), armean_ci(jj,2,:)] = eCI(mustar(:), .95);
        %allmeans(jj,2,:) = mustar;
        
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
        if plotAcorr
            figure(acorrFigh); hold on;
            [acorr, lags] = xcorr(ts(epochi), corrWidth, 'unbiased');
            plot(lags, acorr, 'Color', [.5 .5 .5]);
            med_vel(jj,2) = median(abs(diff(ts(epochi))));
        end
    else
        slopes(jj,2) = NaN; slope_ci(jj,2) = NaN;
    end
    
    
    if ~isempty(occr_trials{jj})
        fprintf('Right Occlusion: ');
        epochi = find(trialCount >= occr_trials{jj}(1) & trialCount <= occr_trials{jj}(end));
        y = ts(epochi);
        fitx = epochi - epochi(1);
        [p, fity, se, ci] = linear_fit_stats(fitx, cum_ts(epochi));
        %plot(epochi, fity,'b', 'LineWidth',1);
        slopes(jj,3) = p(1); slope_ci(jj,3) = ci(1);
        
        % fitting autoregressive (AR with mean) model to timeseries
        [fitp{jj,3}, fitp_ci{jj,3}, e{jj,3}, fit{jj,3}, dwstat(jj,3)] = fitARmodel(ts(epochi), AR_width);
        param_sum = sum(fitp{jj,3}(2:end));  fprintf('Unit root of fit is %f\n', param_sum);
        ar_means(jj,3) =  fitp{jj,3}(1) ./ (1-sum(fitp{jj,3}(2:end))); %mean of the AR model
        [armean(jj,3), armean_ci(jj,3,:), allmeans(jj,3,:)] = mbBootstrap(nboot, anonARmean, ts(epochi), segL, segSep);
        %[armean(jj,3), armean_ci(jj,3,:)] = mbBootstrap(nboot, @fitCumTsSlope, ts(epochi), segL, segSep);
        % Using stationary bootstrap - just going to get means directly from bootstrapped data
        Bstar = opt_block_length_REV_dec07(ts(epochi)); % gives the optimal length for bootstrap blocks
        [ystar, mustar] = stationaryBootstrap(nboot, ceil(Bstar(1)), y);
        %[armean(jj,3), armean_ci(jj,3,:)] = eCI(mustar(:), .95);
        %allmeans(jj,3,:) = mustar;
        
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
        if plotAcorr
            figure(acorrFigh); hold on;
            [acorr, lags] = xcorr(ts(epochi), corrWidth, 'unbiased');
            plot(lags, acorr, 'b');
            med_vel(jj,3) = median(abs(diff(ts(epochi))));
        end
    else
        slopes(jj,3) = NaN; slope_ci(jj,3) = NaN;
    end
    if ~isempty(occl_trials{jj})
        fprintf('Left Occlusion: ');
        epochi = find(trialCount >= occl_trials{jj}(1) & trialCount <= occl_trials{jj}(end));
        y = ts(epochi);
        fitx = epochi - epochi(1);
        [p, fity, se, ci] = linear_fit_stats(fitx, cum_ts(epochi));
        %plot(epochi, fity,'r', 'LineWidth',1);
        slopes(jj,4) = p(1); slope_ci(jj,4) = ci(1);
        
         % fitting autoregressive (AR with mean) model to timeseries
        [fitp{jj,4}, fitp_ci{jj,4}, e{jj,4}, fit{jj,4}, dwstat(jj,4)] = fitARmodel(ts(epochi), AR_width);
        param_sum = sum(fitp{jj,4}(2:end));  fprintf('Unit root of fit is %f\n', param_sum);
        ar_means(jj,4) =  fitp{jj,4}(1) ./ (1-sum(fitp{jj,4}(2:end))); %mean of the AR model
        [armean(jj,4), armean_ci(jj,4,:), allmeans(jj,4,:)] = mbBootstrap(nboot, anonARmean, ts(epochi), segL, segSep);
        %[armean(jj,4), armean_ci(jj,4,:)] = mbBootstrap(nboot, @fitCumTsSlope, ts(epochi), segL, segSep);
        % Using stationary bootstrap - just going to get means directly from bootstrapped data
        Bstar = opt_block_length_REV_dec07(ts(epochi)); % gives the optimal length for bootstrap blocks
        [ystar, mustar] = stationaryBootstrap(nboot, ceil(Bstar(1)), y);
        %[armean(jj,4), armean_ci(jj,4,:)] = eCI(mustar(:), .95);
        %allmeans(jj,4,:) = mustar;
        
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
        if plotAcorr
            figure(acorrFigh); hold on;
            [acorr, lags] = xcorr(ts(epochi), corrWidth, 'unbiased');
            plot(lags, acorr, 'r');
            med_vel(jj,4) = median(abs(diff(ts(epochi))));
        end
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
ar_means_p = ar_means(:,[4 1 2 3]);
% These are using the newer of the fitting routines
armean_ci_p = armean_ci(:,[4 1 2 3],:);
armean_p = armean(:,[4 1 2 3]);
allmeans_p = allmeans(:, [4 1 2 3], :);

%% Plots the fitted biases across conditions

%Simple significance
sig = zeros(size(armean));
armean_sd = NaN * zeros(size(armean_p));
for ii = 1:size(armean,1)
    ctl = squeeze(cat(3, allmeans_p(ii,2,:), allmeans_p(ii,3,:)));
    ctl = ctl(~isnan(ctl));
    armean_p(ii,2) = nanmean(ctl);
    armean_p(ii,3) = nanmean(ctl);
    armean_sd(ii,2) = nanstd(ctl);
    occ = squeeze(allmeans_p(ii,1,:));
    occ = occ(~isnan(occ));
    armean_sd(ii,1) = nanstd(occ);
    if ~isempty(occ) && ranksum(ctl, occ)    
%     if (armean_ci_p(ii,1,2) < armean_ci_p(ii,3,1))
        sig(ii,1) = 1;
    end
    occ = squeeze(allmeans_p(ii,4,:));
    occ = occ(~isnan(occ));
    armean_sd(ii,4) = nanstd(occ);
    if ~isempty(occ) && ranksum(ctl, occ)    
    %if (armean_ci_p(ii,4,1) > armean_ci_p(ii,3,2))
        sig(ii,4) = 1;
    end
end
figure; ah = axes; hold on;
x = [1, 2.25, 2.75, 4];
x = repmat(x, size(ar_means,1), 1);
plotConnectedCategoricalPoints(ah, x', armean_p', sig');
%addErrorBarsAsym(ah, x', armean_p', armean_ci_p(:,:,1)', armean_ci_p(:,:,2)', 'k', .05); 
addErrorBars(ah, x', armean_p', armean_sd', 'k', .05); 




