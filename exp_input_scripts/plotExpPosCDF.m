% Take a structure loaded in the workspace (perMouseData) and come up with a CDF of nose position around
% the trail.

thresh_dist = 40; %px for now
vel_thresh = 40; %mm/sec
mm_conv = .862; 

% all_dists = [];
% for ii = 1:length(exp.vids)
%     segs = exp.vids(ii).getFollowingSegments(1,thresh_dist);
%     odist = [];
%     for jj = 1:size(segs,1)
%         inds = segs(jj,1):segs(jj,2);
%         temp = exp.vids(ii).orthogonalDistFromTrail(inds,1);
%         odist = 
%     end 
% end



%% Plot the following biases of the mice, horizontal position relative to the trail

% We need to define the trials for which the mice are occluded. Since they are occluded in different
% orders this needs to be specified on a per mouse basis.
%ctl_trials = {42:72, 46:76, 43:73, 44:74, 26:45, 26:45, 26:45, 26:45};
%occr_trials = {73:103, 77:107, 74:111, 75:104, 26:35, 26:36, 26:35, 26:35}; %each cell is a mouse
%occl_trials = {16:25, 16:25, 16:25, 16:25, 16:25, 16:25, 16:25, 16:25};
%ctl2_trails = {100:103
mm_conv = .862; %mm per px
nMice = length(perMouseData);
nRows = ceil(sqrt(nMice));
fh = figure; 
sah = axes; hold on;
anal_trials = [5:20];
params = NaN*zeros(nMice, 2);
for ii = 1:nMice
    %sah = subplot(nRows, nRows, ii); hold on;
    rew_free = perMouseData(ii).rew_dists_from_trail(ctl_trials{ii});
    %rew_free = perMouseData(ii).rew_dists_from_trail(anal_trials);
    if(~isempty(rew_free))
        params(ii,:) = plotDistFromTrailDistribution(rew_free, thresh_dist, '', sah);
        
    end
end