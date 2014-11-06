% Let's get the frames when the animal sniffed
for ii=1:length(exp.resp)
    frame_t = exp.camTrig(ii).time(exp.camTrig(ii).frameInds);
    sniff_times = exp.resp(ii).time(exp.resp(ii).sniffVect);
    sniffFrames = zeros(length(sniff_times),1)*NaN;
    for jj=1:length(sniff_times)
        tempf =  find(frame_t >= sniff_times(jj), 1, 'first');
        if ~isempty(tempf)
            sniffFrames(jj) = tempf;
        else 
            sniffFrames(jj) = length(frame_t);
        end
    end
    exp.resp(ii).sniff_times = sniff_times;
    exp.resp(ii).sniff_frames = sniffFrames;
    exp.resp(ii).sniff_pos = exp.vids(ii).nosePos(sniffFrames,:);
end

%%
% Let's plot  
vidi = 1:length(exp.vids);
for ii = 1:length(vidi)
    %exp.vids(ii).plotNosePosition([]);
    exp.vids(ii).plotNosePosition([]);
    %vids(ii).plotFollowing([],15,'');
    %exp.vids(ii).plotPosition([exp.resp(ii).sniffFrames(exp.resp(ii).vidSniffs)], ah, 0, 'b', '.');
end

%% Reset the fcArea
vidi = 1:length(exp.vids);
for ii = 1:length(vidi)
    exp.vids(ii).setFrameCountArea();
end

%% fix nose positions at 0
vidi = 1:length(exp.vids);
for ii = 1:length(vidi)
    np = exp.vids(ii).nosePos;
    %corner = find(np(:,1) <= 1 | np(:,1) >= 1280 | np(:, 2) <= 1 | np(:,2) >= 1024);
    corner = find(np(:,1) >= 1260 | np(:,1) <= 2); 
    exp.vids(ii).nosePos(corner, :) = NaN;
    exp.vids(ii).computeVelocity([]);
    exp.vids(ii).save;
end

%%
for ii = 1:length(exp.vids)
    disp(exp.vids(ii).videoFN);
    exp.vids(ii).blobID(100,:)
    %exp.vids(ii).fcPeriod = 60;
    exp.vids(ii).save;
    
end

%%
for ii = 1:length(exp.vids)
    exp.vids(ii).fcPeriod = 50;
    exp.vids(ii).save;
end

%% Check for missing frames
for ii = 1:length(exp.vids)
    exp.vids(ii).isMissingFrame();
end

%% Plot the velocities of the body and nose
figure; 
for ii = 1:length(exp.vids)
    subplot (3, ceil(length(exp.vids)/3), ii);
    bc = exp.vids(ii).bodyCOM;
    np = exp.vids(ii).nosePos;
    lh = plot(bc(:,1), bc(:,2)); 
    hold on; lh2 = plot(np(:,1), np(:,2), 'g'); 
end

figure; 
for ii = 1:length(exp.vids)
    subplot (3, ceil(length(exp.vids)/3), ii);
    %np = exp.vids(ii).findNose(1:exp.vids(ii).nFrames);
    np = exp.vids(ii).nosePos;
    nf = exp.vids(ii).nFrames;
    disp('\n');
    disp(['Number of frames without nose: ' num2str(sum(isnan(np(:,1))))]);
    exp.vids(ii).computeVelocity(1:exp.vids(ii).nFrames);
    bv = sqrt(sum(exp.vids(ii).bodyVel .^2, 2));
    nv = sqrt(sum(exp.vids(ii).noseVel.^2, 2));
    disp(['Number of frames without nose velocity: ' num2str(sum(isnan(nv)))]);
    %lh = plot(exp.vids(ii).times(:), bv, exp.vids(ii).times(:), nv); 
    lh = plot(exp.vids(ii).times(:), nv, 'Color', 'g', 'LineWidth', 1); 
    %set(lh(2), 'LineWidth', 1);%plot(exp.vids(ii).times(:), bv);
    disp(['Total frames: ' num2str(nf)]);
    disp(['Percent without nose: ' num2str(sum(isnan(np(:,1)))./nf*100)]);
    disp(['Percent without nose velocity: ' num2str(sum(isnan(nv)./nf*100))]);
    %exp.vids(ii).fcPeriod = 60;
    %exp.vids(ii).isMissingFrame();
    %exp.vids(ii).plotNosePosition([]);
end

%% Look at the correlation between velocity (body and nose) with sniff frequency
mm_conv = 1.16; %mm/px linear
f1= figure;
f2 = figure;
f3 = figure;
all_freq = []; all_noseVel_filt = []; all_noseVel = [];
mov_thresh = 50;
for ii = 1:length(exp.resp)
    % select the frames to analyze - when the animal is following the trail
    sniffFrames = exp.resp(ii).sniffFrames(exp.resp(ii).vidSniffs);
    sniffDists = exp.vids(ii).orthogonalDistFromTrail(sniffFrames, 1);
    followi = sniffDists <= 20;
    sniffFrames = sniffFrames(followi); 
    noseVel = exp.vids(ii).noseVel; %get the nose/body velocities
    bodyVel = exp.vids(ii).bodyVel(:,1);
    noseVel_filt = gaussianFilter(noseVel, 3, 'conv'); %smoother versions - vels tend to look messy
    bodyVel_filt = gaussianFilter(bodyVel, 3, 'conv'); 
    noseVel = noseVel(sniffFrames); noseVel_filt = noseVel_filt(sniffFrames);%select the frames
    bodyVel = bodyVel(sniffFrames); bodyVel_filt = bodyVel_filt(sniffFrames);
    sniffFrames = exp.camTrig(ii).frameRange(sniffFrames); %index frames collected rather than those analysed
    freq = exp.resp(ii).sniffFreq(exp.camTrig(ii).frameInds(sniffFrames));
    %noseVel_filt = noseVel_filt ./ nanmax(noseVel_filt); %normalize
    %bodyVel_filt = bodyVel_filt ./ nanmax(bodyVel_filt); 
    norm_freq = freq ./ nanmax(freq);
    % unit conversions - from px/frame to mm/sec
    fps = exp.vids(ii).frameRate;
    noseVel = mm_conv * fps * noseVel; noseVel_filt = mm_conv * fps * noseVel_filt; 
    bodyVel = mm_conv * fps * bodyVel; bodyVel_filt = mm_conv * fps * bodyVel_filt; 
    % Let's try to see what things look like excluding the really low vel points
    moving = noseVel_filt >= mov_thresh;
    freq = freq(moving); noseVel = noseVel(moving); bodyVel = bodyVel(moving);
    noseVel_filt = noseVel_filt(moving);
    % save them
    all_freq = cat(1, all_freq, freq(:));
    all_noseVel_filt = cat(1,all_noseVel_filt, noseVel_filt(:));
    all_noseVel = cat(1,all_noseVel, noseVel(:));
    %correlation/fit
    %[p, S] = polyfit(noseVel, freq, 1);
    %fitx = [0 max(noseVel)];
    %fity = polyval(p, fitx);
    %r = corr(noseVel, freq);
    figure(f1);
    subplot (3, ceil(length(exp.vids)/3), ii);
    plot(noseVel_filt, freq, 'b.', 'MarkerSize', 12);
    %plot(noseVel, freq, 'g.', 'MarkerSize', 12);
    %hold on; plot(bodyVel, freq,'mx');
    %text(50, 8, ['y = ' num2str(p(1)) 'x + ' num2str(p(2)) ', r = ' num2str(r)]); 
    %plot(fitx, fity, 'b', 'LineWidth', 1);
    xlim([50 400]); ylim([5 18]);
    xlabel('Vel (mm/sec)'); ylabel('Sniff freq (Hz)');
    lg_str = {'Nose Vel', 'Body Vel'};
    %legend(lg_str, 'Location', 'SouthEast');
    
    %figure(f2);
    %subplot(3, ceil(length(exp.vids)/3), ii);
    %plot(1:15, 1:15,'k--'); hold on;
    %plot(bodyVel, noseVel, 'r.'); xlabel('Body Vel'); ylabel('Nose Vel');
    %xlim([0 400]); ylim([0 400]);
     
    figure(f3);
    subplot(3, ceil(length(exp.vids)/3), ii);
    x = 1:length(freq);
    %plot(x, noseVel_filt, 'b-', x, bodyVel_filt, 'm-', x, 50*norm_freq, 'r');
    lh = plot(x, noseVel_filt, 'b-', x, freq*3, 'r');
    set(lh, 'LineWidth', 1);
    xlabel('Sniff Number during following'); ylabel('Velocity(mm/s) and Freq (Hz * 3)');
end

f4 = figure;
plot(all_noseVel_filt, all_freq, 'b.', 'MarkerSize', 12);
xlim([50 400]); ylim([5 18]);
xlabel('Vel (mm/sec)'); ylabel('Sniff freq (Hz)');
%correlation/fit
[p, S] = polyfit(all_noseVel_filt, all_freq, 1)
fitx = [mov_thresh max(all_noseVel_filt)];
fity = polyval(p, fitx);
r = corr(all_noseVel, all_freq);
hold on; plot(fitx, fity, 'r');
% Build a binned curve
bs = 25;
vel_bin = mov_thresh:bs:400;
binned_freq = zeros(length(vel_bin)-1,1);
for jj = 1:(length(vel_bin)-1)
    bi = all_noseVel_filt >= vel_bin(jj) & all_noseVel_filt < vel_bin(jj+1);
    binned_freq(jj) = nanmean(all_freq(bi));
end
vel_x = vel_bin(1:end-1) + bs/2;
hold on; plot(vel_x, binned_freq, 'r', 'LineWidth',1);

 %% Plotting respiration frequency as a color over the positions
ft = 1:length(exp.resp);
%ft = 1:3
%make a fake sniff hist to get a good, yet constant color scale
fake_sniff = 13 + 3*randn(1000,1);
fake_sniff = cat(1, fake_sniff, linspace(0,20,40)');
[cm, ~, cvals] = getIndexedColors('jet', fake_sniff,1);
cvals(end) = 50;
%ft = 4;
%fr = 1:1500;

for ii = ft
    fr = 1:exp.vids(ii).nFrames;
    exp.vids(ii).plotPosition(fr, [], 0, 'k', '.');
    pos = exp.resp(ii).sniffPos;
    sniffFrames = exp.resp(ii).sniffFrames(exp.resp(ii).vidSniffs);
    [sniffFrames, sfi,~] = intersect(sniffFrames, fr);
    sniffFrames = exp.camTrig(ii).frameRange(sniffFrames);
    freq = exp.resp(ii).sniffFreq(exp.camTrig(ii).frameInds(sniffFrames)); 
    pos = pos(sfi, :);
    %minfreq = min(freq)
    %maxfreq = max(freq)
    for jj=1:length(freq)
        cind = find(cvals >= freq(jj),1,'first');
        if ~isempty(cind)
            line('Parent',gca,'Xdata',pos(jj,1),'Ydata',pos(jj,2),'Marker','.','MarkerSize',8,'Color',cm(cind,:));
        end
    end
    % to get the scale
    figure; colormap(cm);
    pcolor([freq, freq]);
    colorbar;
end


%%
for ii = 1:length(exp.vids)
    %exp.vids(ii).detectPaths(1,1,1);
    %exp.vids(ii).refinePaths(1);
    %exp.vids(ii).refinePaths(2);
    exp.vids(ii).save;
end

%%
disp('Checking for missing frames in videos');
for ii = 1:length(exp.vids)
    missing = exp.vids(ii).isMissingFrame();
    if missing
        disp(sprintf('%s is missing a frame.', exp.vids(ii).videoFN));
    end
end
%% Let's try to look at the following and sniffing together.

  ft = [3 4 6];
  %ft = 7;
  %fr = 330:900;
  for ii = ft
     fr = 1:exp.vids(ii).nFrames;
     exp.vids(ii).plotFollowing(fr, 15,'');
     fr = intersect(fr, exp.resp(ii).sniffFrames(exp.resp(ii).vidSniffs));
     exp.vids(ii).plotNosePosition(fr, gca, 0, 'y', '.');
  end
  
 %%
 nose_mean = NaN*ones(size(nose_px,1),1);
 for ii = 1:size(nose_px,1)
    nose_mean(ii) = nanmean(nose_px(ii,:));
 end
 
 %% Plot position/velocity for videos
 fh = figure;
 mm_conv = 1.16;
 for ii =1:length(exp.vids)
     exp.vids(ii).plotVelocity([],'nose', 2);
     
     figure(fh);
     subplot (3, ceil(length(exp.vids)/3), ii);
     nv = exp.vids(ii).noseVel;
     nv_filt = gaussianFilter(nv, 3, 'conv');
     nv_filt = nv_filt * mm_conv * exp.vids(ii).frameRate;
     t = exp.vids(ii).times;
     plot(t, nv_filt, 'b', 'LineWidth', 1);
     xlabel('Time (sec)');
     ylabel('Nose Velocity');
     %ylim([0 300]);
 end
 
 %%
 for ii =1:length(exp.vids)
     exp.vids(ii).makePathsSkel;
 end
 %%
 for ii =1:length(exp.vids)
     exp.vids(ii).plotOrientation([]);
 end
 
 %%
 traj = []; dir = []; wind = [-15:60];
 for ii =1:length(exp.vids)
    exp.vids(ii).makePathsSkel();
    [traj_tmp, dir_tmp] = noseTrajectories(exp.vids(ii));
    if ii == 1
        traj = traj_tmp;
        dir = dir_tmp;
    else
        traj = cat(1,traj, traj_tmp);
        dir = cat(1, dir, dir_tmp);
    end
 end
 
 neg = (dir==-1);
 mean_neg_traj = nanmean(traj(neg,:));
 pos = (dir==1);
 mean_pos_traj = nanmean(traj(pos,:));
 
 figure;
 plot(wind, traj(neg,:), 'b'); hold on;
 plot(wind, mean_neg_traj, 'b', 'LineWidth',2);
 plot(wind, traj(pos,:), 'm');
 plot(wind, mean_pos_traj, 'm','LineWidth', 2);
 
 %% 
 for ii=1:length(exp.resp)
    figure;
    sniff = exp.resp(ii).value;
    sniff = sniff-mean(sniff);
    sniff = sniff./max(sniff);
    plot(exp.resp(ii).time, sniff*2+20, 'k');
    hold on;
    plot(exp.resp(ii).time, exp.resp(ii).sniffFreq, 'r');
 end
 
    
     
     
 