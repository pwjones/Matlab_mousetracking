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
vidi = 1:length(vids);
%vidi = 1:2;
%vidi = 24:28;
for ii = 1:length(vidi)
    %exp.vids(ii).plotNosePosition([]);
    %vids(ii).plotPosition([]);
    vids(ii).plotFollowing([],15,'');
    %exp.vids(ii).plotPosition([exp.resp(ii).sniffFrames(exp.resp(ii).vidSniffs)], ah, 0, 'b', '.');
end
%%
for ii = 1:length(exp.vids)
    disp(exp.vids(ii).videoFN);
    exp.vids(ii).blobID(100,:)
    %exp.vids(ii).fcPeriod = 60;
    exp.vids(ii).save;
    
end


%% Plot the velocities of the body and nose
figure; 
for ii = 1:length(exp.vids)
    subplot (3, ceil(length(exp.vids)/3), ii);
    bc = exp.vids(ii).bodyCOM;
    lh = plot(bc(:,1), bc(:,2)); 
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
    lh = plot(exp.vids(ii).times(:), bv, exp.vids(ii).times(:), nv); 
    set(lh(2), 'LineWidth', 1);%plot(exp.vids(ii).times(:), bv);
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
    noseVel_filt = noseVel_filt ./ nanmax(noseVel_filt); %normalize
    bodyVel_filt = bodyVel_filt ./ nanmax(bodyVel_filt); 
    norm_freq = freq ./ nanmax(freq);
    % unit conversions - from px/frame to mm/sec
    fps = exp.vids(ii).frameRate;
    noseVel = mm_conv * fps * noseVel; noseVel_filt = mm_conv * fps * noseVel_filt; 
    bodyVel = mm_conv * fps * bodyVel; bodyVel_filt = mm_conv * fps * bodyVel_filt; 
    % Let's try to see what things look like excluding the really low vel points
    %moving = noseVel >= 75;
    %freq = freq(moving); noseVel = noseVel(moving); bodyVel = bodyVel(moving);
    %correlation/fit
    [p, S] = polyfit(noseVel, freq, 1);
    fitx = [0 max(noseVel)];
    fity = polyval(p, fitx);
    r = corr(noseVel, freq);
    figure(f1);
    subplot (3, ceil(length(exp.vids)/3), ii);
    plot(noseVel_filt, freq, 'bo');
    %hold on; plot(bodyVel, freq,'mx');
    %text(50, 8, ['y = ' num2str(p(1)) 'x + ' num2str(p(2)) ', r = ' num2str(r)]); 
    %plot(fitx, fity, 'b', 'LineWidth', 1);
    %xlim([0 400]); ylim([5 Inf]);
    xlabel('Vel (mm/sec)'); ylabel('Sniff freq (Hz)');
    lg_str = {'Nose Vel', 'Body Vel'};
    legend(lg_str, 'Location', 'SouthEast');
    
    figure(f2);
    subplot(3, ceil(length(exp.vids)/3), ii);
    plot(1:15, 1:15,'k--'); hold on;
    plot(bodyVel, noseVel, 'r.'); xlabel('Body Vel'); ylabel('Nose Vel');
    xlim([0 400]); ylim([0 400]);
    
    figure(f3);
    subplot(3, ceil(length(exp.vids)/3), ii);
    x = 1:length(norm_freq);
    %plot(x, noseVel_filt, 'b-', x, bodyVel_filt, 'm-', x, 50*norm_freq, 'r');
    plot(x, noseVel_filt, 'b-', x, 50*norm_freq, 'r');
    xlabel('Sniff Number during following'); ylabel('Velocity(mm/s) and Freq (norm. 0-50)');
end

%% Plotting respiration frequency as a color over the positions
ft = 1:length(exp.resp);
%ft = 3;
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
    [cm cinds] = getIndexedColors('jet', freq, 1);
    %minfreq = min(freq)
    %maxfreq = max(freq)
    for jj=1:length(freq)
        line('Parent',gca,'Xdata',pos(jj,1),'Ydata',pos(jj,2),'Marker','.','MarkerSize',8,'Color',cm(cinds(jj),:));
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
 for ii =1:length(exp.vids)
     exp.vids(ii).plotVelocity([],'nose', 0);
     
     figure(fh);
     subplot (3, ceil(length(exp.vids)/3), ii);
     nv = exp.vids(ii).noseVel;
     t = exp.vids(ii).times;
     plot(t, nv, 'b', 'LineWidth', 1);
     xlabel('Time (sec)');
     ylabel('Nose Velocity');
 end
 
 
 