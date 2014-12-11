% Take a structure loaded in the workspace and come up with a CDF of nose position around the trail

thresh_dist = 20; %px for now
vel_thresh = 40; %mm/sec
mm_conv = .862; 

all_dists = [];
for ii = 1:length(exp.vids)
    segs = exp.vids(ii).getFollowingSegments(1,thresh_dist);
    odist = [];
    for jj = 1:size(segs,1)
        inds = segs(jj,1):segs(jj,2);
        temp = exp.vids(ii).orthogonalDistFromTrail(inds,1);
        odist = 
end