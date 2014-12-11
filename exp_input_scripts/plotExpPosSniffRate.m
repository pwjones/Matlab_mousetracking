% Plotting respiration frequency as a color over the positions

ft = 1:length(exp.resp);
%ft = 1:3
ft=3;
%make a fake sniff hist to get a good, yet constant color scale
%fake_sniff = 13 + 3*randn(1000,1);
fake_sniff = (8:.05:17)';
fake_sniff = cat(1, fake_sniff, linspace(0,20,20)');
%[cm, cm_inds, cvals] = getIndexedColors('jet', fake_sniff,1);
cvals(end) = 50;
cm_range = [10 18];

for ii = ft
    fr = 1:exp.vids(ii).nFrames;
    exp.vids(ii).plotPosition(fr, [], 0, 'k', '.');
    pos = exp.resp(ii).sniffPos;
    sniffFrames = exp.resp(ii).sniffFrames(exp.resp(ii).vidSniffs);
    [sniffFrames, sfi,~] = intersect(sniffFrames, fr);
    sniffFrames_abs = exp.camTrig(ii).frameRange(sniffFrames);
    freq = exp.resp(ii).sniffFreq(exp.camTrig(ii).frameInds(sniffFrames_abs)); 
    pos = exp.vids(ii).nosePos(sniffFrames, :);
    
    [cm, cm_inds, cvals] = getIndexedColors('jet',freq, 0, cm_range);
    %minfreq = min(freq)
    %maxfreq = max(freq)
    for jj=1:length(freq)
        %cind = find(cvals >= freq(jj),1,'first');
        line('Parent',gca,'Xdata',pos(jj,1),'Ydata',pos(jj,2),'Marker','.','MarkerSize',8,'Color',cm(cm_inds(jj),:));
    end
    % to get the scale
    figure; colormap(cm);
    pcolor([cm_range(1):cm_range(2); cm_range(1):cm_range(2)]);
    colorbar;
end