% popTurningAnal.m

% population analysis of the turning behavior
mm_conv = .862; %mm per px
nMice = length(perMouseData);
nRows = ceil(sqrt(nMice));
fh = figure; hold on;
dl = [-1, 1];
wind = -15:15;
selx = 16:46;
centerx = 31;

overlay_fig = figure; ah = axes; hold on;
figure; % The figure for all of the trajectory plots
for ii = 1:nMice
    % free the data from the structure
    rew_free = perMouseData(ii).turning_traj(ctl_trials{ii});
    dirs = perMouseData(ii).turning_dir(ctl_trials{ii});
    subplot(nRows, ceil(nMice/nRows), ii);
    xlim([min(wind) max(wind)]);
    set(gca, 'YDir', 'reverse');
    
    free = [];
    for kk = 1:length(rew_free)
        free = cat(2, free, rew_free{kk});
    end
    
    if (ii==1) 
        %turning distances from trail by trial, dir, epoch, mouse
        turn_dists = NaN*zeros(size(free,2), 2, 4, nMice); 
    end
    dir_arr = cell2mat(dirs);
    mean_turns = []; se_turns = [];
    for jj = 1:length(dl)
        sel = (dir_arr == dl(jj))';
        nTurn = sum(sel);
        mean_turns(:,jj) = nanmean(free(selx,sel),2);
        turn_dists(1:nTurn, jj, 1, ii) = free(centerx, sel);  
        if ~isempty(mean_turns)
            ci = bootci(1000, {@nanmean, free(selx,sel)'},'type','cper');
            plot_err_poly_asym(gca, wind', mean_turns(:,jj), ci', 'k', [.5 .5 .5], 1);
            hold on;
            plot(ah, wind', mean_turns(:,jj), 'k-');
        end
    end
    
    % NEED TO INCLUDE THE OTHER NON OCCLUDED PERIOD
    
    % Second Unoccluded Period 
    rew_free2 = perMouseData(ii).turning_traj(ctl2_trials{ii});
    dirs = perMouseData(ii).turning_dir(ctl2_trials{ii});
    
    free2 = [];
    for kk = 1:length(rew_free2)
        free2 = cat(2, free2, rew_free2{kk});
    end
    
    dir_arr = cell2mat(dirs);
    mean_turns = [];
    for jj = 1:length(dl)
        if ~isempty(free2)
            sel = (dir_arr == dl(jj))';
            nTurn = sum(sel);
            mean_turns(:,jj) = nanmean(free2(selx,sel),2);
            turn_dists(1:nTurn, jj, 2, ii) = free2(centerx, sel);  
            ci = bootci(1000, {@nanmean, free2(selx,sel)'},'type','cper');
            plot_err_poly_asym(gca, wind', mean_turns(:,jj), ci', 'k', [.5 .5 .5], 1);
        end
    end
    
    % Right Occlusions 
    rew_occr = perMouseData(ii).turning_traj(occr_trials{ii});
    dirs = perMouseData(ii).turning_dir(occr_trials{ii});
    
    occr = [];
    for kk = 1:length(rew_occr)
        occr = cat(2, occr, rew_occr{kk});
    end
    
    dir_arr = cell2mat(dirs);
    mean_turns = [];
    for jj = 1:length(dl)
        if ~isempty(occr)
            sel = (dir_arr == dl(jj))';
            nTurn = sum(sel);
            mean_turns(:,jj) = nanmean(occr(selx,sel),2);
            turn_dists(1:nTurn, jj, 3, ii) = occr(centerx, sel);  
            ci = bootci(1000, {@nanmean, occr(selx,sel)'},'type','cper');
            plot_err_poly_asym(gca, wind', mean_turns(:,jj), ci', 'b', [.717 .855 1], 1);
        end
    end
    
    % Left occlusions
    rew_occl = perMouseData(ii).turning_traj(occl_trials{ii});
    dirs = perMouseData(ii).turning_dir(occl_trials{ii});
    
    occl = [];
    for kk = 1:length(rew_occl)
        occl = cat(2, occl, rew_occl{kk});
    end
    
    dir_arr = cell2mat(dirs);
    mean_turns = [];
    for jj = 1:length(dl)
        if ~isempty(occl)
            sel = (dir_arr == dl(jj))';
            nTurn = sum(sel);
            mean_turns(:,jj) = nanmean(occl(selx,sel),2);
            turn_dists(1:nTurn, jj, 4, ii) = occl(centerx, sel);  
            ci = bootci(1000,{@nanmean, occl(selx,sel)'},'type','cper');
            plot_err_poly_asym(gca, wind', mean_turns(:,jj), ci', 'r', [1 .78 .847], 1);
        end
    end
    
end

zi = turn_dists == 0; 
turn_dists(zi) = NaN;

%% Plotting the mean shifts and computing significance
figure; axes;
plot([1.8,4.2], [0 0], 'k--');
xlim([1.8, 4.2]);
hold on;
set(gca, 'Ydir', 'reverse');

pb = 0;
histx = -20:.5:20;
mean_turn_dists = squeeze(nanmean(turn_dists));
h = NaN*zeros(size(mean_turn_dists));
sig_cond = zeros(size(mean_turn_dists));
mean_turn_dists_p = NaN * zeros(size(mean_turn_dists));
for ii = 2:size(turn_dists,3)
    for jj = 1:2
        for kk = 1:nMice
            x1 = turn_dists(:,jj,1,kk);
            x1 = x1(~isnan(x1));
            x2 = turn_dists(:,jj,ii,kk);
            x2 = x2(~isnan(x2));
            if pb
               figure;
               y1 = histc(x1, histx);
               y2 = histc(x2, histx);
               bar(histx', [y1, y2]);
            end
            
            if ~isempty(x1) && ~isempty(x2)
                [p,h] = ranksum(x1,x2);
                sig_cond(jj,ii,kk) = h;
                mean_turn_dists_p(jj,ii,kk) = mean_turn_dists(jj,ii,kk) - mean_turn_dists(jj,1,kk);
%                 (if (h == 1)
%                     plot(ii, mean_turn_dists(jj,ii,kk) - mean_turn_dists(jj,1,kk), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 12);
%                 else
%                     plot(ii, mean_turn_dists(jj,ii,kk) - mean_turn_dists(jj,1,kk), 'ko', 'MarkerSize', 12);
%                 end
%                 hold on;
                
            end
        end
    end
end

mean_turn_dists_p = mean_turn_dists_p(:,2:end, :);
sig_cond = sig_cond(:,2:end,:);
%reorder the dims
mean_turn_dists_p = permute(mean_turn_dists_p, [2 3 1]);
sig_cond = permute(sig_cond, [2 3 1]);
% collapse across turn directions
sig_cond = reshape(sig_cond, size(sig_cond,1), []);
mean_turn_dists_p = reshape(mean_turn_dists_p, size(mean_turn_dists_p,1), []);
% reorder the conditions
sig_cond = sig_cond([3 1 2], :);
mean_turn_dists_p = mean_turn_dists_p([3 1 2], :);
%now plot it
figure; ah = axes; hold on;
x = [1.5, 2, 2.5];
plot([1.4, 2.6], [0 0], 'k--');
xlim([1.4, 2.6]);
hold on;
set(gca, 'Ydir', 'reverse', 'TickDir', 'Out');

x = repmat(x, size(mean_turn_dists_p,2), 1);
plotConnectedCategoricalPoints(ah, x', mean_turn_dists_p, sig_cond);

