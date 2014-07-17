%% Plot the time on trail over time.
%plot(1:length(fake_prop), fake_prop *100,'Color', [.25 .25 .25]); hold on;
figure; ah1 = axes; hold on;
figure; ah2 = axes; hold on;
title('Trail Fraction Followed per Trial');
nMice = length(perMouseData)/nConc;
nRows = ceil(sqrt(nMice));
filtn = 3; boxcar = ones(1,filtn)./filtn; %define the averaging filter kernel
for ii = 1:(nMice*nConc)
    %subplot(nRows, nRows, ii); %square, many panels
    nTrials = length(perMouseData(ii).rew_dists);
    plot(ah1, perMouseData(ii).rew_prop*100, 'Color', gcolor{genotype(ii)}); hold on;
    plot(ah1, perMouseData(ii).dist_prop*100, 'r'); hold on;
    plot(ah2, (perMouseData(ii).rew_prop/(perMouseData(ii).rew_prop+perMouseData(ii).dist_prop))*100, 'Color', gcolor{genotype(ii)}); hold on;
    samps_cut = floor(length(boxcar)/2);
    vi = (1+samps_cut):(length(perMouseData(ii).rew_prop)-samps_cut)
    xl = length(vi)+samps_cut;
    %plot([1 xl], 100*[mean(perMouseData(ii).rew_prop) mean(perMouseData(ii).rew_prop)], 'k--');
    %plot([1 xl], [35 35], '--k'); hold on;
    %fake_prop_filt = conv(fake_prop, boxcar,'valid');
    rew_prop_filt = conv(perMouseData(ii).rew_prop, boxcar,'valid');
    dist_prop_filt = conv(perMouseData(ii).dist_prop, boxcar,'valid');
    %vi = (1+samps_cut):(samps_cut+length(rew_prop_filt));
    %flh = plot(vi, fake_prop_filt*100,'LineWidth',2, 'Color', [.25 .25 .25]);
    plot(ah1, vi, rew_prop_filt*100, 'Color', gcolor{genotype(ii)}, 'LineWidth',3);
    plot(ah1, vi, dist_prop_filt*100, 'r', 'LineWidth',3);
    set(gca, 'TickDir','out', 'fontsize', 16);
    xlabel('Trial #','FontSize', 18);
    ylabel('% Time on Trail','FontSize', 18);
    title(mouse_names{ii});
    % Let's do a slightly different measure - ratio of following rewarded/total
    plot(ah2, vi(:), (rew_prop_filt./(rew_prop_filt+dist_prop_filt))*100, 'Color', gcolor{genotype(ii)}, 'LineWidth',3);
    set(gca, 'TickDir','out', 'fontsize', 16);
    xlabel('Trial #','FontSize', 18);
    ylabel('% Time Following Rewarded Trail','FontSize', 18);
    title(mouse_names{ii});
    %legend({'Rewarded Trail', 'Distracter Trail'});
    %title('Proportion of Time Following Trails');
end
