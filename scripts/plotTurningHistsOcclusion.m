%function all_mags_dd = plotTurningHists(sniffData)
% %function all_mags_dd = plotTurningHists(sniffData)
%
% perMouseData(1)  = 
% rew_prop: [178�1 double]
% dist_prop: [178�1 double]
% fake_prop: [178�1 double]
% rew_dists: {178�1 cell}
% distract_dists: {178�1 cell}
% total_frames: [178�1 double]
% frame_rate: [178�1 double]
% rew_dists_from_trail: {178�1 cell}
% distract_dists_from_trail: {178�1 cell}
% total_turning: [178�1 double]
% rew_dists_from_trail_persect: {178�1 cell}
% distract_dists_from_trail_persect: {178�1 cell}
% rew_trail_area: [178�1 double]
% dist_trail_area: [178�1 double]
% rew_propFollowed: [178�1 double]
% dist_propFollowed: [178�1 double]
% dir_nums: [178�1 double]
% nose_trajectories: {178�1 cell}
% traj_window: [1�61 double]
% traj_dir: {178�1 cell}
% turning_traj: {178�1 cell}
% turning_dir: {178�1 cell} % dir - direction of the turn. -1: rightward, 1: leftward
% file_names: {178�1 cell}
% nose_vel: {178�1 cell}
% body_vel: {178�1 cell}

%% Plot the overall turning directions
fallCohortList;

nMice = length(perMouseData);
nRows = ceil(sqrt(nMice));
conds = {'ctl', 'occr', 'occl'}; %control (open), right occluded, left occluded
allDirs = cell(length(conds),1); allTurnPos = cell(length(conds),1); 
allTraj = cell(length(conds),1); allDirsTowards = cell(length(conds),1);
turni = 31; %middle of the window where the mouse turns (poor coding to hard code this - I'm ashamed)
for cc = 1:length(conds)
    cond = eval([conds{cc} '_trials']); %variable specifying the trials for each condition
    for ii = 1:nMice
        vids = cond{ii}; %the trials that correspond to the current condition
        for jj = vids
            temp = perMouseData(ii).turning_dir{jj};
            allDirs{cc} = cat(1, allDirs{cc}, temp);
            %turni = find(perMouseData(ii).traj_window == 0, 1, 'first');
            traj = perMouseData(ii).turning_traj{jj};
            allTraj{cc} = cat(2, allTraj{cc}, traj);
            allTurnPos{cc} = cat(1, allTurnPos{cc}, traj(turni,:)');
            turnTowards = temp .* traj(turni,:)';
            turnTowards = turnTowards ./ abs(turnTowards);
            allDirsTowards{cc} = cat(1, allDirsTowards{cc}, turnTowards);
        end
    end
end

%% Now plot the histograms of turning direction as a function of distance
% from the trail.
edges = -20:2:20;
figure; ah = axes; 
plotstyles = {'k-', 'r-', 'b-'};
turnTowardsBinned = cell(length(conds),1);
for cc = 1:length(conds)
    [ah, turnTowardsBinned{cc}] = makeTurningHistogram(allDirsTowards{cc}, allTurnPos{cc}, edges, ah, plotstyles{cc});
    hold on;
    turnTowards = allDirsTowards{cc};
    overallTurnTowards(cc) = sum(turnTowards(turnTowards == 1)) ./ length(turnTowards);
end
legend(conds);



%% This is the code that actually makes/plots the histograms
% turnp = NaN*zeros(numel(edges),1);
% 
% [N,bin] = histc(allTurnPos, edges);
% nbins = length(N);
% propTurnsLeft = zeros(nbins, 1);
% %turnTowardsProp = zeros(1,nbins);
% for ii = 1:nbins
%     temp = allDirs(bin == ii); % turns at each distance
%     % Direction is coded as follows, -1: rightward, 1: leftward
%     propTurnsLeft(ii) = sum(temp == 1) ./ length(temp);
%     %turnTowards(ii) = sum(temp(temp == 1)); %The number of turns towards the trail
%     %turnTowardsProp(ii, jj-1) = turnTowards(ii)/numel(temp(~isnan(temp)));
%     %temp = turnsLeft(bin == ii);
%     %nTurnsLeft = sum(temp > 0);
%     %pTurnsLeft(ii, jj-1) = nTurnsLeft./numel(temp(~isnan(temp)));
% end
% 
% figure;
% stairs(edges, propTurnsLeft);
% xlabel('Distance from Trail (mm)');
% ylabel('Prop Turns Leftwards');
% %legend({'Sniff -2', 'Sniff -1','Sniff 0'});
% xlim([min(edges) max(edges)]);



%% population analysis of the turning behavior
% mm_conv = .862; %mm per px
% nMice = length(perMouseData);
% nRows = ceil(sqrt(nMice));
% fh = figure; hold on;
% dl = [-1, 1];
% wind = -15:15;
% selx = 16:46;
% centerx = 31;
% 
% overlay_fig = figure; ah = axes; hold on;
% figure; % The figure for all of the trajectory plots
% for ii = 1:nMice
%     % free the data from the structure
%     rew_free = perMouseData(ii).turning_traj(ctl_trials{ii});
%     dirs = perMouseData(ii).turning_dir(ctl_trials{ii});
%     subplot(nRows, ceil(nMice/nRows), ii);
%     xlim([min(wind) max(wind)]);
%     set(gca, 'YDir', 'reverse');
%     
%     free = [];
%     for kk = 1:length(rew_free)
%         free = cat(2, free, rew_free{kk});
%     end
%     
%     if (ii==1) 
%         %turning distances from trail by trial, dir, epoch, mouse
%         turn_dists = NaN*zeros(size(free,2), 2, 4, nMice); 
%     end
%     dir_arr = cell2mat(dirs);
%     mean_turns = []; se_turns = [];
%     for jj = 1:length(dl)
%         sel = (dir_arr == dl(jj))';
%         nTurn = sum(sel);
%         mean_turns(:,jj) = nanmean(free(selx,sel),2);
%         turn_dists(1:nTurn, jj, 1, ii) = free(centerx, sel);  
%         if ~isempty(mean_turns)
%             ci = bootci(1000, {@nanmean, free(selx,sel)'},'type','cper');
%             plot_err_poly_asym(gca, wind', mean_turns(:,jj), ci', 'k', [.5 .5 .5], 1);
%             hold on;
%             plot(ah, wind', mean_turns(:,jj), 'k-');
%         end
%     end
%     
%     % NEED TO INCLUDE THE OTHER NON OCCLUDED PERIOD
%     
%     % Second Unoccluded Period 
%     rew_free2 = perMouseData(ii).turning_traj(ctl2_trials{ii});
%     dirs = perMouseData(ii).turning_dir(ctl2_trials{ii});
%     
%     free2 = [];
%     for kk = 1:length(rew_free2)
%         free2 = cat(2, free2, rew_free2{kk});
%     end
%     
%     dir_arr = cell2mat(dirs);
%     mean_turns = [];
%     for jj = 1:length(dl)
%         if ~isempty(free2)
%             sel = (dir_arr == dl(jj))';
%             nTurn = sum(sel);
%             mean_turns(:,jj) = nanmean(free2(selx,sel),2);
%             turn_dists(1:nTurn, jj, 2, ii) = free2(centerx, sel);  
%             ci = bootci(1000, {@nanmean, free2(selx,sel)'},'type','cper');
%             plot_err_poly_asym(gca, wind', mean_turns(:,jj), ci', 'k', [.5 .5 .5], 1);
%         end
%     end
%     
%     % Right Occlusions 
%     rew_occr = perMouseData(ii).turning_traj(occr_trials{ii});
%     dirs = perMouseData(ii).turning_dir(occr_trials{ii});
%     
%     occr = [];
%     for kk = 1:length(rew_occr)
%         occr = cat(2, occr, rew_occr{kk});
%     end
%     
%     dir_arr = cell2mat(dirs);
%     mean_turns = [];
%     for jj = 1:length(dl)
%         if ~isempty(occr)
%             sel = (dir_arr == dl(jj))';
%             nTurn = sum(sel);
%             mean_turns(:,jj) = nanmean(occr(selx,sel),2);
%             turn_dists(1:nTurn, jj, 3, ii) = occr(centerx, sel);  
%             ci = bootci(1000, {@nanmean, occr(selx,sel)'},'type','cper');
%             plot_err_poly_asym(gca, wind', mean_turns(:,jj), ci', 'b', [.717 .855 1], 1);
%         end
%     end
%     
%     % Left occlusions
%     rew_occl = perMouseData(ii).turning_traj(occl_trials{ii});
%     dirs = perMouseData(ii).turning_dir(occl_trials{ii});
%     
%     occl = [];
%     for kk = 1:length(rew_occl)
%         occl = cat(2, occl, rew_occl{kk});
%     end
%     
%     dir_arr = cell2mat(dirs);
%     mean_turns = [];
%     for jj = 1:length(dl)
%         if ~isempty(occl)
%             sel = (dir_arr == dl(jj))';
%             nTurn = sum(sel);
%             mean_turns(:,jj) = nanmean(occl(selx,sel),2);
%             turn_dists(1:nTurn, jj, 4, ii) = occl(centerx, sel);  
%             ci = bootci(1000,{@nanmean, occl(selx,sel)'},'type','cper');
%             plot_err_poly_asym(gca, wind', mean_turns(:,jj), ci', 'r', [1 .78 .847], 1);
%         end
%     end
%     
% end
% 
% zi = turn_dists == 0; 
% turn_dists(zi) = NaN;

%%
% 
% % Figure 1 - Distance changes, turning probabilities
% % ---------------------------------------------------
% turnDirections = sniffData.turnDirections./abs(sniffData.turnDirections); %convert to signed 1's.
% edges = -9.5:.5:9.5;
% for jj = 1:size(sniffData.dist_diffs,2)
%     [N,bin] = histc(sniffData.dist_diffs(:,jj), edges);
%     nbins = length(N);
%     turn_counts = zeros(1,nbins);
%     turnTowards= zeros(1,nbins);
%     turnTowardsProp = zeros(1,nbins);
%     for ii = 1:nbins
%         turn_counts(ii) = sum(sniffData.bturn(bin == ii));
%         temp = turnDirections(bin == ii); % turns after each concentration change
%         turnTowards(ii) = sum(temp(temp == 1)); %The number of turns towards the trail
%         turnTowardsProp(ii) = turnTowards(ii)/turn_counts(ii);
%         %headings = mouseHeading(bin == ii);
%         %meanHeadings(ii) = mean(headings);
%         %stdHeadings(ii) = std(headings);
%     end
%     turnp(:,jj) = turn_counts(:)./N(:);
% end
% turnp(isnan(turnp)) = 0; %replace possible NaNs
% % similar, but doing it using the turning triggered data.
% for ii = 1:nbins
%     for jj = 1:3
%         [N_s, bin_s] = histc(sniffData.turnTrig_preTurnDistDiff(:,jj), edges);
%         turns = sniffData.turnTrig_turnDirs(bin_s==ii);
%         turnToP(ii,jj) = sum(turns>0)./numel(turns);
%     end
% end
%  
% % getting the actual sets of delta NDs in order to test the distributions
% % for unimodality
% for jj=1:size(sniffData.dist_diffs, 2)
%    turni = sniffData.bturn == 1;
%    dND(:,jj) = sniffData.dist_diffs(turni,jj); 
% end
% 
% figure;
% subplot(1,3,1);
% %bar(edges(1:end-1), N, 'histc');
% stairs(edges, N);
% hold on;
% stairs(edges, N_s);
% xlabel('\Delta Distance from Trail (mm)');
% ylabel('Counts');
% xlim([-10 10]);
% 
% subplot(1,3,2);
% %bar(edges(1:end-1),turnp, 'histc');
% stairs(edges',turnp);
% xlabel('\Delta Distance from Trail (mm)');
% ylabel('Turning Probability');
% xlim([-10 10]);
% ylim([0 .5]);
% 
% subplot(1,3,3);
% %bar(edges(1:end-1),turnTowardsProp, 'histc');
% %bar(edges, turnTowardsProp, 'histc');
% stairs(edges', turnToP);
% xlabel('\Delta Distance from Trail (mm)');
% ylabel('Prop Turns Towards Trail');
% %legend({'-3 to -4', '-3 to -2','-2 to -1'});
% xlim([-10 10]);
% 
% % Figure 2 - Headings Pre and Post Turn
% % -------------------------------------
% figure;
% subplot(2,3,1);
% cats = [-10 -2.5 2.5 10];
% [Ncat, catbin] = histc(sniffData.turnTrig_preTurnDistDiff(:,3), cats);
% pre = sniffData.turnTrig_preTurnHeadings(:,1);
% post = sniffData.turnTrig_postTurnHeadings(:,1);
% rh = rose2(mod(pre(catbin == 1)+2*pi, 2*pi));
% title('Approach - Pre');
% %---------------
% subplot(2,3,2);
% rh = rose2(mod(pre(catbin == 2)+2*pi, 2*pi));
% title('Small Change - Pre');
% %---------------
% subplot(2,3,3);
% rh = rose2(mod(pre(catbin == 3)+2*pi, 2*pi));
% title('Away - Pre');
% %---------------
% subplot(2,3,4);
% rh = rose2(mod(post(catbin == 1)+2*pi, 2*pi));
% title('Approach - Post Turn');
% %---------------
% subplot(2,3,5);
% rh = rose2(mod(post(catbin == 2)+2*pi, 2*pi));
% title('Small Change - Post Turn');
% %---------------
% subplot(2,3,6);
% rh = rose2(mod(post(catbin == 3)+2*pi, 2*pi));
% title('Away - Post Turn');
% 
% % Figure 3,4 - Turning Magnitude Scatterplot and Histograms
% % -------------------------------------
% plotblue = [.3 .6 1];
% scat_fig = figure;
% % Distance decreases - approaching
% pre_angles = mod(sniffData.turnTrig_preTurnHeadings(catbin == 1)+(2*pi), 2*pi);
% post_angles = mod(sniffData.turnTrig_postTurnHeadings(catbin == 1)+(2*pi), 2*pi);
% plot(pre_angles, post_angles, '^g'); hold on;
% diff_fig = figure; 
% labelah = axes('Position', [ 0 0 1 1], 'Visible', 'off'); 
% text(.3, .1, ' Turning Magnitudes', 'FontSize', 18);
% subplot(2,3,1);
% ang_diff = circ_dist(post_angles, pre_angles);
% rose2(ang_diff);
% title('Intersniff Approach');
% %[~, ang_diff] = rotationDirection(pre_angles, post_angles);
% med_diff = median(abs(ang_diff))
% % Magnitude histograms
% subplot(2,3,4); [N, cent] = hist(abs(ang_diff)); bar(cent, N, 'hist'); 
% xlim([0 pi]); hold on; plot(gca, [med_diff med_diff], [0 max(N) + 10], 'k--');
% text(med_diff, max(N)+3, ['median = ' num2str(med_diff)]);
% set(gca, 'Xtick', [0 pi/2 pi], 'XTickLabel', {'0', 'p/2', 'p'}, 'FontName', 'symbol');
% ha = findobj(gca,'Type','patch'); set(ha, 'FaceColor', plotblue);
% xlabel('Turn Magnitude (radians)', 'FontName', 'Helvetica'); ylabel('Count', 'FontName', 'Helvetica');
% 
% % Distance spanning zero
% pre_angles = mod(sniffData.turnTrig_preTurnHeadings(catbin == 2)+(2*pi), 2*pi);
% post_angles = mod(sniffData.turnTrig_postTurnHeadings(catbin == 2)+(2*pi), 2*pi);
% figure(scat_fig); plot(pre_angles, post_angles, 'xk');
% figure(diff_fig); subplot(2,3,2);
% ang_diff = circ_dist(post_angles, pre_angles);
% rose2(ang_diff);
% title('Intersniff Small Change');
% %[~, ang_diff] = rotationDirection(pre_angles, post_angles);
% med_diff = median(abs(ang_diff))
% % Magnitude histograms
% subplot(2,3,5); [N, cent] = hist(abs(ang_diff)); bar(cent, N, 'hist'); 
% xlim([0 pi]); hold on; plot(gca, [med_diff med_diff], [0 max(N) + 10], 'k--');
% text(med_diff, max(N)+3, ['median = ' num2str(med_diff)]);
% set(gca, 'Xtick', [0 pi/2 pi], 'XTickLabel', {'0', 'p/2', 'p'}, 'FontName', 'symbol');
% ha = findobj(gca,'Type','patch'); set(ha, 'FaceColor', plotblue);
% xlabel('Turn Magnitude (radians)', 'FontName', 'Helvetica'); ylabel('Count', 'FontName', 'Helvetica');
% title('Turns', 'FontName', 'Helvetica');
% % Distance increases - away
% pre_angles = mod(sniffData.turnTrig_preTurnHeadings(catbin == 3)+(2*pi), 2*pi);
% post_angles = mod(sniffData.turnTrig_postTurnHeadings(catbin == 3)+(2*pi), 2*pi);
% figure(scat_fig); plot(pre_angles, post_angles, 'or');
% xlim([0 2*pi]); ylim([0 2*pi]);
% figure(diff_fig); subplot(2,3,3);
% ang_diff = circ_dist(post_angles, pre_angles);
% rose2(ang_diff);
% title('Intersniff Away');
% %[~, ang_diff] = rotationDirection(pre_angles, post_angles);
% med_diff = median(abs(ang_diff))
% % Magnitude histograms
% subplot(2,3,6); [N, cent] = hist(abs(ang_diff)); bar(cent, N, 'hist'); 
% xlim([0 pi]); hold on; plot(gca, [med_diff med_diff], [0 max(N) + 10], 'k--');
% text(med_diff, max(N)+3, ['median = ' num2str(med_diff)]);
% set(gca, 'Xtick', [0 pi/2 pi], 'XTickLabel', {'0', 'p/2', 'p'}, 'FontName', 'symbol');
% ha = findobj(gca,'Type','patch'); set(ha, 'FaceColor', plotblue);
% xlabel('Turn Magnitude (radians)', 'FontName', 'Helvetica'); ylabel('Count', 'FontName', 'Helvetica');
% %axes(labelah);
% %text(.3, .1, 'Turning Magnitudes');
% 
% %% Figure 5 - Dependence of mean turning magnitude on the distance differences and the 
% %             absolute position of the mouse.
% corder = {[.9 .4 .4], [1 0 0], [0 0 0]};
% edges_dd = -10:4:10;
% edges_pos = -14:4:14;
% post = convertMouseHeadingAngles(sniffData.turnTrig_postTurnHeadings, sniffData.turnTrig_postTurnPos, sniffData.turnTrig_postTurnDirToTrail);
% pre = convertMouseHeadingAngles(sniffData.turnTrig_preTurnHeadings, sniffData.turnTrig_preTurnPos, sniffData.turnTrig_preTurnDirToTrail);
% ang_diff = circ_dist(post, pre);
% turn_mag = abs(ang_diff);
% med_mags_dd = NaN*zeros(length(edges_dd),3);
% mean_mags_dd = NaN*zeros(length(edges_dd),3);
% ste_mags_dd = NaN*zeros(length(edges_dd),3);
% med_mags_pos = NaN*zeros(length(edges_pos),3);
% mean_mags_pos = NaN*zeros(length(edges_dd),3);
% ste_mags_pos = NaN*zeros(length(edges_dd),3);
% all_mags_dd = {};
% for jj = 1:3 % different # of sniffs back
%     [N, bin] = histc(sniffData.turnTrig_preTurnDistDiff(:,jj), edges_dd);
%     for ii = 1:length(N)
%         mags = 180/pi*turn_mag(bin==ii);
%         med_mags_dd(ii, jj) = median(mags);
%         mean_mags_dd(ii, jj) = mean(mags);
%         ste_mags_dd(ii, jj) = std(mags)./sqrt(numel(mags));
%         all_mags_dd{ii, jj} = mags;
%     end
%     
%     [N, bin] = histc(sniffData.turnTrig_sniffPos(:,jj+1), edges_pos);
%     for ii = 1:length(N)
%         mags = 180/pi*turn_mag(bin==ii);
%         med_mags_pos(ii, jj) = median(mags);
%         mean_mags_pos(ii, jj) = mean(mags);
%         ste_mags_pos(ii, jj) = std(mags)./sqrt(numel(mags));
%         all_mags_dd{ii, jj} = mags; %uncomment to test position rather than change
%     end
% end
% 
% % Let's test each of these lines for flatness via bootstrap methods
% bootstrap_h = zeros(size(all_mags_dd));
% bootstrap_pthresh = .05 / (size(all_mags_dd,1)-1);
% bootstrap_p = zeros(size(all_mags_dd));
% for jj=1:3 %# of sniffs back
%     comps = [1,2,4,5];
%     for ii=1:length(comps)
%         [bootstrap_h(ii,jj), p] = bootstrapMeanTest(all_mags_dd{3, jj}, all_mags_dd{comps(ii), jj}, bootstrap_pthresh);
%         bootstrap_p(ii,jj) = min(p);
%     end
% end
% 
% % Two way ANOVA on the turning magnitudes
% all_mags_dd = all_mags_dd(1:end-1, :);
% %all_mags_dd = all_mags_dd';
% % buildANOVA2Matrix;
% % disp('Doing an ANOVA on turning magnitudes');
% % [p, table, stats] = anova2(testM, minelem);
% % [c,m,h] = multcompare(stats, 'ctype', 'bonferroni');
% 
% 
% figure;
% subplot(2,2,1);
% plot(edges_dd', med_mags_dd, 'k-o');
% xlabel('\Delta ND (mm)'); ylabel('Median Turn Mag'); xlim([-10 10]);
% subplot(2,2,2); hold on;
% x = edges_dd(1:end-1) + (edges_dd(2)-edges_dd(1))/2;
% for ii = 1:3
%     plot_err_poly(gca, x', mean_mags_dd(1:end-1,ii), ste_mags_dd(1:end-1,ii), corder{ii}, (corder{ii}+[1 1 1])/2, .5);
% end
% xlabel('\Delta ND (mm)'); ylabel('Mean Turn Mag'); xlim([-10 10]);
% subplot(2,2,3);
% plot(edges_pos', med_mags_pos,'k-o');
% xlabel('ND (mm)'); ylabel('Median Turn Mag'); xlim([-15 15]);
% subplot(2,2,4); hold on;
% x = edges_pos(1:end-1) + (edges_pos(2)-edges_pos(1))/2;
% for ii = 1:3
%     plot_err_poly(gca, x', mean_mags_pos(1:end-1,ii), ste_mags_pos(1:end-1,ii), corder{ii}, (corder{ii}+[1 1 1])/2, .5);
% end
% xlabel('ND (mm)'); ylabel('Mean Turn Mag'); xlim([-15 15]);
% %plot(sniffData.turnTrig_preTurnDistDiff(:,3), turn_mag, '.k', 'MarkerSize', 6);
% %subplot(1,3,3);
% %magM = cell2mat_2Dunequal(all_mags{:,3});
% %boxplot(magM, 'notch', 'on');
% 
% %% Figure 6 - Plotting measures of turning relative to sniff location
% % rather than intersniff differences
% 
% turnDirections = sniffData.turnDirections./abs(sniffData.turnDirections); %convert to signed 1's.
% turnPos = sniffData.sniffPos(:,end);
% turnsLeft = turnDirections.*turnPos;
% edges = -15:2:15;
% turnp = NaN*zeros(numel(edges), size(sniffData.sniffPos,2)-1);
% for jj = 2:size(sniffData.sniffPos,2)
%     [N,bin] = histc(sniffData.sniffPos(:,jj), edges);
%     nbins = length(N);
%     turn_counts = zeros(1,nbins);
%     turnTowards= zeros(1,nbins);
%     %turnTowardsProp = zeros(1,nbins);
%     for ii = 1:nbins
%         turn_counts(ii) = sum(sniffData.bturn(bin == ii));
%         temp = turnDirections(bin == ii); % turns after each concentration change
%         turnTowards(ii) = sum(temp(temp == 1)); %The number of turns towards the trail
%         turnTowardsProp(ii, jj-1) = turnTowards(ii)/numel(temp(~isnan(temp)));
%         temp = turnsLeft(bin == ii);
%         nTurnsLeft = sum(temp > 0);
%         pTurnsLeft(ii, jj-1) = nTurnsLeft./numel(temp(~isnan(temp)));
%     end
%     turnp(:,jj-1) = turn_counts(:)./N(:);
% end
% 
% figure;
% subplot(1, 3, 1);
% stairs(edges',turnp);
% xlabel('Distance from Trail (mm)');
% ylabel('Turning Probability');
% xlim([-15 15]);
% ylim([0 .5]);
% 
% subplot(1,3,2);
% stairs(edges', turnTowardsProp);
% xlabel('Distance from Trail (mm)');
% ylabel('Prop Turns Towards Trail');
% %legend({'Sniff -2', 'Sniff -1','Sniff 0'});
% xlim([-15 15]);
% 
% subplot(1,3,3);
% stairs(edges', pTurnsLeft);
% xlabel('Distance from Trail (mm)');
% ylabel('Prop Turns Leftwards');
% %legend({'Sniff -2', 'Sniff -1','Sniff 0'});
% xlim([-15 15]);
% 
% %% Figure 7 - Turning Magnitudes as a function of sniff position
% 
