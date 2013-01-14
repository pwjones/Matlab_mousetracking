distance_comp = cell(2,2); %want to do a reward/distractor, early/late comparison
trial_range = [1:8; 37:44];
for ii = 1:2
    trials = trial_range(ii,:);
    for jj=1:length(trials)
        distance_comp{ii,1} = cat(1, distance_comp{ii,1}, rew_dists{trials(jj)});
        distance_comp{ii,2} = cat(1, distance_comp{ii,2}, distract_dists{trials(jj)});
    end
end
    
counts = cell(2,2); %want to do a reward/distractor, early/late comparison
[counts{1,1}, xbins] = hist(distance_comp{1,1}, 400); counts{1,2} = hist(distance_comp{1,2},xbins);
counts{2,1} = hist(distance_comp{2,1},xbins); counts{2,2} = hist(distance_comp{2,2}, xbins);
figure;
subplot(2,1,1); hold on;
bar(xbins, counts{1,1}, 'g'); hold on; 
bar(xbins, -counts{1,2}, 'r'); hold on;
xlim([0 2000]);
ylim([-500 500]);
title('First 8 Trials');
subplot(2,1,2); hold on;
bar(xbins, counts{2,1}, 'g'); hold on; 
bar(xbins, -counts{2,2}, 'r'); hold on;
xlim([0 2000]);
ylim([-500 500]);
title('Last 8 Trials');
