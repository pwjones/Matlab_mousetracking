nMice = length(perMouseData);
nRows = ceil(sqrt(nMice));
distThresh = 15;

%initialize
AR_width = 12;
fitp = cell(nMice, 4); fitp_ci = cell(nMice, 4); %parameters of model fit
e = cell(nMice, 4); fit = cell(nMice, 4);

figure;
for jj = 1:nMice
    subplot(nRows, nRows, jj); hold on;
    ts = []; trialCount = [];
    mouseDists = perMouseData(jj).rew_dists_from_trail;
    for ii = 1:length(mouseDists)
        % First plot some markers for epoch beginnings
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
        following = mouseDists{ii}(abs(mouseDists{ii}) <= distThresh);
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
        epochi = find(trialCount >= ctl_trials{jj}(1) & trialCount <= ctl_trials{jj}(end));
        %linear trend fit of cumulative
        fitx = epochi - epochi(1);
        [p, fity, se, ci] = linear_fit_stats(fitx, cum_ts(epochi));
        plot(epochi, fity,'Color',[.5 .5 .5], 'LineWidth',1);
        slopes(jj,1) = p(1); slope_ci(jj,1) = ci(1);
        % fitting autoregressive (with mean) model to timeseries
        fprintf('Control: ');
        [fitp{jj,1}, fitp_ci{jj,1}, e{jj,1}, fit{jj,1}] = fitARmodel(ts(epochi), AR_width);
        plot(epochi((AR_width+2):end), cum_ts(epochi(1)+AR_width-1)+cumsum(fit{jj,1}), '--', 'Color',[.5 .5 .5]);
    end
    if ~isempty(ctl2_trials{jj})
        epochi = find(trialCount >= ctl2_trials{jj}(1) & trialCount <= ctl2_trials{jj}(end));
        fitx = epochi - epochi(1);
        [p, fity, se, ci] = linear_fit_stats(fitx, cum_ts(epochi));
        plot(epochi, fity,'Color',[.5 .5 .5], 'LineWidth',1);
        slopes(jj,2) = p(1); slope_ci(jj,2) = ci(1);
        % fitting autoregressive (with mean) model to timeseries
        fprintf('Control2: ');
        [fitp{jj,2}, fitp_ci{jj,2}, e{jj,2}, fit{jj,2}] = fitARmodel(ts(epochi), AR_width);
        plot(epochi((AR_width+2):end), cum_ts(epochi(1)+AR_width-1)+cumsum(fit{jj,2}), 'k--');
    else
        slopes(jj,2) = NaN; slope_ci(jj,2) = NaN;
    end
    if ~isempty(occr_trials{jj})
        epochi = find(trialCount >= occr_trials{jj}(1) & trialCount <= occr_trials{jj}(end));
        fitx = epochi - epochi(1);
        %[p, S] = polyfit(fitx, cum_ts(epochi),1);
        %fity = polyval(p,fitx,S);
        [p, fity, se, ci] = linear_fit_stats(fitx, cum_ts(epochi));
        plot(epochi, fity,'b', 'LineWidth',1);
        slopes(jj,3) = p(1); slope_ci(jj,3) = ci(1);
        % fitting autoregressive (with mean) model to timeseries
        fprintf('Right Occlusion: ');
        [fitp{jj,3}, fitp_ci{jj,3}, e{jj,3}, fit{jj,3}] = fitARmodel(ts(epochi), AR_width);
        plot(epochi((AR_width+2):end), cum_ts(epochi(1)+AR_width-1)+cumsum(fit{jj,3}), 'k--');
    else
        slopes(jj,3) = NaN; slope_ci(jj,3) = NaN;
    end
    if ~isempty(occl_trials{jj})
        epochi = find(trialCount >= occl_trials{jj}(1) & trialCount <= occl_trials{jj}(end));
        fitx = epochi - epochi(1);
        %[p, S] = polyfit(fitx, cum_ts(epochi),1);
        %fity = polyval(p,fitx,S);
        [p, fity, se, ci] = linear_fit_stats(fitx, cum_ts(epochi));
        plot(epochi, fity,'r', 'LineWidth',1);
        slopes(jj,4) = p(1); slope_ci(jj,4) = ci(1);
        % fitting autoregressive (with mean) model to timeseries
        fprintf('Left Occlusion: ');
        [fitp{jj,4}, fitp_ci{jj,4}, e{jj,4}, fit{jj,4}] = fitARmodel(ts(epochi), AR_width);
        plot(epochi((AR_width+2):end), cum_ts(epochi(1)+AR_width-1)+cumsum(fit{jj,4}), 'k--');
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
for ii = 1:numel(fitp)
    temp = fitp{ii};
    tci = fitp_ci{ii};
    if ~isempty(temp)
        ar_slopes(ii) = temp(1);
        slope_lb(ii) = tci(1,1);
        slope_ub(ii) = tci(1,2);
    end
end
