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

%% Plot the following proportions for the trails  
vidi = 1:length(exp.vids);
figure; ah = axes; hold on;
for ii = 1:length(vidi)
    exp.vids(ii).plotFollowingTimecourse(20, ah);
end
ylim([0 1]);

%%
vidi = 1:length(exp.vids(ii));
for ii = 1:length(vidi)    
    exp.vids(ii).computeVelocity([]);
end

%% Reset the fcArea
vidi = 1:length(exp.vids);
for ii = 1:length(vidi)
    exp.vids(ii).setFrameCountArea();
end

%% delete nose positions around edges
vidi = 1:length(exp.vids);
for ii = 1:length(vidi)
    np = NaN*ones(exp.vids(ii).nFrames, 2);
    nb = exp.vids(ii).noseblob;
    for jj = 1:exp.vids(ii).nFrames
        if ~isnan(nb(jj))
            np(jj, :) = exp.vids(ii).areas(jj,nb(jj)).Centroid;
        end
    end
    edge = find(np(:,1) <= 15 | np(:,1) >= 1265 | np(:, 2) <= 15 | np(:,2) >= 1010)
    %corner = find(np(:,1) >= 1260 | np(:,1) <= 2); 
    exp.vids(ii).nosePos(edge, :) = NaN;
    exp.vids(ii).noseblob(edge) = NaN;
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
    exp.vids(ii).save;
end


%%
for ii = 1:length(exp.vids)
    exp.vids(ii).fcPeriod = 50;
    exp.vids(ii).save;
end

%% Check for missing frames
for ii = 1:length(exp.vids)
    exp.vids(ii).fcArea = [1221, 959, 1280, 999];
    %exp.vids(ii).fcArea
    exp.vids(ii).isMissingFrame();
end

%% change the conversions
for ii = 1:length(exp.vids)
   exp.vids(ii).mm_conv = .862;
   exp.vids(ii).save;
    
end

%% change the conversions
for ii = 1:length(exp.vids)
   exp.vids(ii).findFollowingTurns([], 1, 20, -20:20);
end

%% Plot the velocities of the body and nose

vel_con = .862 * exp.vids(1).frameRate; %mm/px linear
vel_fig = figure;
pos_fig = figure; 
for ii = 1:length(exp.vids)
    subplot (3, ceil(length(exp.vids)/3), ii);
    bc = exp.vids(ii).bodyCOM;
    np = exp.vids(ii).nosePos;
    lh = plot(bc(:,1), bc(:,2)); 
    hold on; lh2 = plot(np(:,1), np(:,2), 'Color', [0 .7 0], 'LineWidth', 1); 
end

pb = 0;
for ii = 1:length(exp.vids)
    figure(vel_fig);
    subplot (3, ceil(length(exp.vids)/3), ii);
    %np = exp.vids(ii).findNose(1:exp.vids(ii).nFrames);
    np = exp.vids(ii).nosePos;
    nf = exp.vids(ii).nFrames;
    disp('\n');
    disp(['Number of frames without nose: ' num2str(sum(isnan(np(:,1))))]);
    exp.vids(ii).computeVelocity(1:exp.vids(ii).nFrames);
    bv = sqrt(sum(exp.vids(ii).bodyVel .^2, 2)) * vel_con;
    nv = sqrt(sum(exp.vids(ii).noseVel.^2, 2)) * vel_con;
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
    toohigh = find(nv > 1000);
    if ~isempty(toohigh) && pb
        figure(pos_fig);
        subplot (3, ceil(length(exp.vids)/3), ii);
        plot(exp.vids(ii).nosePos(toohigh,1), exp.vids(ii).nosePos(toohigh,2), 'rx');
    end
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

%% Plot the sniffing and the detected
ft = 1
for ii = ft
    figure; 
    hold on; 
    plot(exp.resp(ii).vidTime, exp.resp(ii).value_filt, 'k');
    sniffi = find(exp.resp(ii).sniffVect > .5);
    plot(exp.resp(ii).vidTime(sniffi), exp.resp(ii).value_filt(sniffi), 'ro');
end
    
%% Plot the velocity and the sniff rate with each other.
ft = 1
mm_conv = .862;
for ii = ft
    v_scale = mm_conv*exp.vids(ii).frameRate;
    figure;
    nv = exp.vids(ii).noseVel;
    mov_t = exp.vids(ii).times;
    nv_filt = gaussianFilter(nv, 1.5, 'conv');
    
    [AX,H1,H2] = plotyy(exp.resp(ft).vidTime, exp.resp(ft).sniffFreq, mov_t, nv_filt*v_scale, 'plot');
    axes(AX(1)); hold on;
    axes(AX(2)); hold on;
    plot(AX(2), [mov_t(1) mov_t(end)], [50 50], 'g--'); %threshold for including in analysis
    rh = plotSpikeRasters(AX(1), exp.resp(ft).vidTime, exp.resp(ft).sniffVect', 17);
    %set(rh, 'Color', 'k', 'LineWidth', 1);
    set(get(AX(1),'Ylabel'),'String','Sniff Rate (Hz)') ;
    set(get(AX(2),'Ylabel'),'String','Nose Velocity (mm/sec)');
    xlabel('Time (sec)') 
    set(AX(1),'Xlim', [mov_t(1) mov_t(end)]);
    set(AX(2), 'Xlim', [mov_t(1) mov_t(end)]);
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
 mm_conv = .862;
 for ii =1:length(exp.vids)
     exp.vids(ii).plotVelocity([],'nose', 2);
     
     figure(fh);
     subplot (3, ceil(length(exp.vids)/3), ii);
     nv = exp.vids(ii).noseVel;
     nv_filt = gaussianFilter(nv, 1, 'conv');
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
 
    
     
     
 