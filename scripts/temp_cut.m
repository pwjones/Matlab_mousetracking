%% Plot the time on trail over time.
%plot(1:length(fake_prop), fake_prop *100,'Color', [.25 .25 .25]); hold on;
% figure; ah1 = axes; hold on;
% figure; ah2 = axes; hold on;
% title('Trail Fraction Followed per Trial');
% nMice = length(perMouseData)/nConc;
% nRows = ceil(sqrt(nMice));
% filtn = 3; boxcar = ones(1,filtn)./filtn; %define the averaging filter kernel
% for ii = 1:(nMice*nConc)
%     %subplot(nRows, nRows, ii); %square, many panels
%     nTrials = length(perMouseData(ii).rew_dists);
%     plot(ah1, perMouseData(ii).rew_prop*100, 'Color', gcolor{genotype(ii)}); hold on;
%     plot(ah1, perMouseData(ii).dist_prop*100, 'r'); hold on;
%     plot(ah2, (perMouseData(ii).rew_prop/(perMouseData(ii).rew_prop+perMouseData(ii).dist_prop))*100, 'Color', gcolor{genotype(ii)}); hold on;
%     samps_cut = floor(length(boxcar)/2);
%     vi = (1+samps_cut):(length(perMouseData(ii).rew_prop)-samps_cut)
%     xl = length(vi)+samps_cut;
%     %plot([1 xl], 100*[mean(perMouseData(ii).rew_prop) mean(perMouseData(ii).rew_prop)], 'k--');
%     %plot([1 xl], [35 35], '--k'); hold on;
%     %fake_prop_filt = conv(fake_prop, boxcar,'valid');
%     rew_prop_filt = conv(perMouseData(ii).rew_prop, boxcar,'valid');
%     dist_prop_filt = conv(perMouseData(ii).dist_prop, boxcar,'valid');
%     %vi = (1+samps_cut):(samps_cut+length(rew_prop_filt));
%     %flh = plot(vi, fake_prop_filt*100,'LineWidth',2, 'Color', [.25 .25 .25]);
%     plot(ah1, vi, rew_prop_filt*100, 'Color', gcolor{genotype(ii)}, 'LineWidth',3);
%     plot(ah1, vi, dist_prop_filt*100, 'r', 'LineWidth',3);
%     set(gca, 'TickDir','out', 'fontsize', 16);
%     xlabel('Trial #','FontSize', 18);
%     ylabel('% Time on Trail','FontSize', 18);
%     title(mouse_names{ii});
%     % Let's do a slightly different measure - ratio of following rewarded/total
%     plot(ah2, vi(:), (rew_prop_filt./(rew_prop_filt+dist_prop_filt))*100, 'Color', gcolor{genotype(ii)}, 'LineWidth',3);
%     set(gca, 'TickDir','out', 'fontsize', 16);
%     xlabel('Trial #','FontSize', 18);
%     ylabel('% Time Following Rewarded Trail','FontSize', 18);
%     title(mouse_names{ii});
%     %legend({'Rewarded Trail', 'Distracter Trail'});
%     %title('Proportion of Time Following Trails');
% end


% figure; hold on;
% for jj = 1:length(perMouseData)
%     med_dist_filt = [conv( perMouseData(jj).med_dist(:,1), boxcar(:),'valid') conv( perMouseData(jj).med_dist(:,2), boxcar(:),'valid')]; 
%     plot(perMouseData(jj).med_dist(:,1), '-', 'LineWidth', .5,'Color', gcolor{genotype(jj)}); hold on;
%     plot(perMouseData(jj).med_dist(:,2), '-', 'Color', 'r', 'LineWidth', .5);
%     vi = (1+samps_cut):(samps_cut+size(med_dist_filt,1));
%     plot(vi, med_dist_filt(:,1),'LineWidth', 2,'Color', gcolor{genotype(jj)}); 
%     plot(vi, med_dist_filt(:,2), 'r', 'LineWidth', 2);    
% end
% legend({'Rewarded Trail', 'Distracter Trail'});
% xlabel('Trial #', 'FontSize', 14); ylabel('Median Following Distance (px)', 'FontSize', 14);

%% Plot the time on trail over time.
% %plot(1:length(fake_prop), fake_prop *100,'Color', [.25 .25 .25]); hold on;
% figure; ah1 = axes; hold on;
% figure; ah2 = axes; hold on;
% figure; ah3 = axes; hold on;
% %title('Trail Fraction Followed per Trial');
% nMice = length(perMouseData)/nConc;
% nRows = ceil(sqrt(nMice));
% filtn = 3; boxcar = ones(1,filtn)./filtn; %define the averaging filter kernel
% for ii = 1:(nMice*nConc)
%     %subplot(nRows, nRows, ii); %square, many panels
%     nTrials = length(perMouseData(ii).rew_dists);
%     plot(ah1, perMouseData(ii).rew_prop*100, 'Color', gcolor{genotype(ii)}); hold on;
%     plot(ah1, perMouseData(ii).dist_prop*100, 'r'); hold on;
%     propFollowing = perMouseData(ii).rew_prop./(perMouseData(ii).rew_prop+perMouseData(ii).dist_prop);
%     plot(ah2, propFollowing*100, 'Color', gcolor{genotype(ii)}); hold on;
%     samps_cut = floor(length(boxcar)/2);
%     vi = (1+samps_cut):(length(perMouseData(ii).rew_prop)-samps_cut)
%     xl = length(vi)+samps_cut;
%     %plot([1 xl], 100*[mean(perMouseData(ii).rew_prop) mean(perMouseData(ii).rew_prop)], 'k--');
%     %plot([1 xl], [35 35], '--k'); hold on;
%     %fake_prop_filt = conv(fake_prop, boxcar,'valid');
%     rew_prop_filt = conv(perMouseData(ii).rew_prop, boxcar,'valid');
%     dist_prop_filt = conv(perMouseData(ii).dist_prop, boxcar,'valid');
%     %vi = (1+samps_cut):(samps_cut+length(rew_prop_filt));
%     %flh = plot(vi, fake_prop_filt*100,'LineWidth',2, 'Color', [.25 .25 .25]);
%     plot(ah1, vi, rew_prop_filt*100, 'Color', gcolor{genotype(ii)}, 'LineWidth',3);
%     plot(ah1, vi, dist_prop_filt*100, 'r', 'LineWidth',3);
%     set(gca, 'TickDir','out', 'fontsize', 16);
%     xlabel('Trial #','FontSize', 18);
%     ylabel('% Time on Trail','FontSize', 18);
%     title(mouse_names{ii});
%     % Let's do a slightly different measure - ratio of following rewarded/total
%     plot(ah2, vi(:), (rew_prop_filt./(rew_prop_filt+dist_prop_filt))*100, 'Color', gcolor{genotype(ii)}, 'LineWidth',3);
%     set(gca, 'TickDir','out', 'fontsize', 16);
%     xlabel('Trial #','FontSize', 18);
%     ylabel('% Time Following Rewarded Trail','FontSize', 18);
%     title(mouse_names{ii});
%     plot(ah3, ones(length(propFollowing),1)*ii, propFollowing, 'LineStyle','.', 'Color', gcolor{genotype(ii)}, 'LineWidth',1);
%     plot(ah3, ii, nanmean(propFollowing), 'LineStyle', 'x', 'Color', gcolor{genotype(ii)}, 'MarkerSize', 12,'LineWidth',2);
%     set(gca, 'TickDir','out', 'fontsize', 16);
%     xlabel('Mouse','FontSize', 18);
%     ylabel('% Time Following Rewarded Trail','FontSize', 18);
%     title(mouse_names{ii});
%     
% end


% %% Plot the time on trail over time.
% figure; ah1 = axes; hold on;
% xlabel('Trial #','FontSize', 18);
% ylabel('% Time on Trail','FontSize', 18);
% figure; ah2 = axes; hold on;
% xlabel('Trial #','FontSize', 18);
% ylabel('% Time Following Rewarded Trail','FontSize', 18);
% figure; ah3 = axes; hold on;
% xlabel('Concentration','FontSize', 18);
% ylabel('Prop Following Rewarded','FontSize', 18);
% nMice = length(perMouseData)/nConc;
% nRows = ceil(sqrt(nMice));
% filtn = 3; boxcar = ones(1,filtn)./filtn; %define the averaging filter kernel
% mean_prop_follow = []; std_prop_follow = []; n_prop_follow = [];
% cl = {[0 0 1], [0 1 0], [0 0 0]}; offset = [-.05 .05];
% rew_colors = blackGradColormap([0 1 0], nConc+1); rew_colors = rew_colors(2:end,:);
% dist_colors = blackGradColormap([1 0 0], nConc+1); dist_colors = dist_colors(2:end,:);
% for jj = 1:nMice
%     for kk=1:nConc
%         ii= (jj-1)*nConc + (kk-1) + 1;
%         %subplot(nRows, nRows, ii); %square, many panels
%         nTrials = length(perMouseData(ii).rew_dists);
%         plot(ah1, perMouseData(ii).rew_prop*100, 'Color', rew_colors(kk,:)); hold on;
%         plot(ah1, perMouseData(ii).dist_prop*100, 'Color', dist_colors(kk,:)); hold on;
%         % Let's also do a slightly different measure - ratio of following rewarded/total
%         propTimeFollowingRew = perMouseData(ii).rew_prop./(perMouseData(ii).rew_prop+perMouseData(ii).dist_prop);
%         mean_prop_follow(kk,jj) = nanmean(propTimeFollowingRew); std_prop_follow(kk,jj) = nanstd(propTimeFollowingRew);
%         n_prop_follow(kk,jj) = nTrials;
%         plot(ah2, propTimeFollowingRew*100, 'Color', rew_colors(kk,:)); hold on;
%         plot(ah3, (kk+offset(jj))*ones(nTrials,1), propTimeFollowingRew, '.', 'Color', cl{jj});
%         samps_cut = floor(length(boxcar)/2);
%         vi = (1+samps_cut):(length(perMouseData(ii).rew_prop)-samps_cut)
%         xl = length(vi)+samps_cut;
%         %plot([1 xl], 100*[mean(perMouseData(ii).rew_prop) mean(perMouseData(ii).rew_prop)], 'k--');
%         %plot([1 xl], [35 35], '--k'); hold on;
%         %fake_prop_filt = conv(fake_prop, boxcar,'valid');
%         rew_prop_filt = conv(perMouseData(ii).rew_prop, boxcar,'valid');
%         dist_prop_filt = conv(perMouseData(ii).dist_prop, boxcar,'valid');
%         %vi = (1+samps_cut):(samps_cut+length(rew_prop_filt));
%         %flh = plot(vi, fake_prop_filt*100,'LineWidth',2, 'Color', [.25 .25 .25]);
%         plot(ah1, vi, rew_prop_filt*100, 'Color', rew_colors(kk,:), 'LineWidth',3);
%         plot(ah1, vi, dist_prop_filt*100, 'Color', dist_colors(kk,:), 'LineWidth',3);
%         set(ah1, 'TickDir','out', 'fontsize', 16);
%         plot(ah2, vi(:), (rew_prop_filt./(rew_prop_filt+dist_prop_filt))*100, 'Color', rew_colors(kk,:), 'LineWidth',3);
%         set(ah2, 'TickDir','out', 'fontsize', 16);
%     end
% end
% % Let's produce a plot of the proportion of time on trails spent following the rewarded trail, as a function of concentration
% 
% for jj=1:nMice
%     plot(ah3, mean_prop_follow(:,jj),'o-', 'Color', cl{jj}); hold on;
%     addErrorBars(ah3, 1:nConc, mean_prop_follow(:,jj), std_prop_follow(:,jj)./sqrt(n_prop_follow(:,jj)), cl{jj}, .1);
% en


%% 
for ii = 1:length(vids)
    np = vids(ii).nosePos;
    corner = find(np(:,1) >= 1206 & np(:, 2) >= 954);
    vids(ii).nosePos(corner, :) = NaN;
    vids(ii).computeVelocity([]);
    vids(ii).save;
end

%%

% %% Plot the time on trail over time - plot by day
% figure; hold on;
% title('Trail Fraction Followed per Day');
% nMice = length(perMouseData);
% nRows = ceil(sqrt(nMice));
% max_days = length(folders);
% for ii = 1:nMice
%     nTrails = length(folder_nums{ii}); % added 10/14 - want to plot by day rather than by trial
%     days = sort(unique(folder_nums{ii})); %the day labels
%     rew_prop = []; dist_prop = [];
%     for jj = 1:length(days)
%         sel = folder_nums{ii} == days(jj); 
%         rew_prop(jj) = nanmean(perMouseData(ii).rew_prop(sel)) * 100;
%         dist_prop(jj) = nanmean(perMouseData(ii).dist_prop(sel)) * 100;
%     end
%     plot(days, rew_prop, 'g-', 'LineWidth', 1); hold on;
%     plot(days, dist_prop, 'r-', 'LineWidth', 1);
%     set(gca, 'TickDir','out', 'fontsize', 16);
%     xlabel('Training Day','FontSize', 18);
%     ylabel('% Time on Trail','FontSize', 18);
%     title(mouse_names{ii});
%     
%     all_days{ii} = days;
%     all_rew_prop{ii} = rew_prop;
%     all_dist_prop{ii} = dist_prop;
% end
% % Now we need to add an overall mean line
% mean_rew_prop = zeros(max_days, 1)*NaN; mean_dist_prop = zeros(max_days, 1)*NaN;
% for ii = 1:max_days
%    rew_prop = zeros(nMice,1) * NaN; dist_prop = zeros(nMice,1) * NaN;
%    for jj = 1:nMice
%       fi = find(all_days{jj} == ii); % find the matching day index for this mouse
%       rew_prop(jj) = all_rew_prop{jj}(fi); % select the rewarded prop corresponding to the day
%       dist_prop(jj) = all_dist_prop{jj}(fi); 
%    end
%    mean_rew_prop(ii) = nanmean(rew_prop);
%    mean_dist_prop(ii) = nanmean(dist_prop);
% end
% plot(1:max_days, mean_rew_prop, 'g-', 'LineWidth', 2);
% plot(1:max_days, mean_dist_prop, 'r-', 'LineWidth', 2);