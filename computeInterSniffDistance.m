function [dists, vels] = computeInterSniffDistance(exp)
% function dists = computeInterSniffDistance(exp)

mm_conv = .862;
ntrial = length(exp.vids);
dists = []; vels=[];

for ii = 1:ntrial
    diffSniffPos = diff(exp.resp(ii).sniffPos);
    sniffDists = [NaN; sqrt(sum(diffSniffPos.^2, 2))] * mm_conv; %convert to mm from px
    sniffFrames = exp.resp(ii).sniffFrames(exp.resp(ii).vidSniffs);
    sniffNoseVels = exp.vids(ii).noseVel(sniffFrames) * mm_conv * exp.vids(ii).frameRate;
    
    dists = cat(1, dists, sniffDists);
    vels = cat(1, vels, sniffNoseVels);
end