% Just go through the tracked videos and plot the nose positions

vidi = 1:length(vids);
%vidi = 1:2;
for ii = 1:length(vidi)
    %exp.vids(ii).plotNosePosition([]);
    vids(ii).plotPosition([]);
end