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
%     if (armean_ci_p(ii,1,2) < armean_ci_p(ii,3,1))
    if ~isempty(occ)
        [p,h] = ranksum(ctl, occ);    
        if h
            sig(ii,1) = 1;
        end
    end
    occ = squeeze(allmeans_p(ii,4,:));
    occ = occ(~isnan(occ));
    armean_sd(ii,4) = nanstd(occ);
    %if (armean_ci_p(ii,4,1) > armean_ci_p(ii,3,2))
    if ~isempty(occ) && ranksum(ctl, occ)   
        [p,h] = ranksum(ctl, occ);    
        if h
            sig(ii,4) = 1;
        end
    end
end
figure; ah = axes; hold on;
x = [1, 2.25, 2.75, 4];
x = repmat(x, size(ar_means,1), 1);
plotConnectedCategoricalPoints(ah, x', armean_p', sig');
%addErrorBarsAsym(ah, x', armean_p', armean_ci_p(:,:,1)', armean_ci_p(:,:,2)', 'k', .05); 
addErrorBars(ah, x', armean_p', armean_sd', 'k', .05); 