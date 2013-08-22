function ret = loadMouseTracking(list_file, file_range, following_thresh)
% function loadMouseTracking(list_file, file_range)
% 
% Function to load a set of processed videos, their saved MouseTracker
% objects. The set is defined in the LIST_FILE, of which subsets can be
% selected using the FILE_RANGE argument. 

base_fname = '/Users/pwjones/Movies/mouse_training/';
%following_thresh = 20; %px, the distance from the trail the animal can get before it's counted as not following
class_func = @MouseTrackerKF;
%%%%%%%%%%%%%%%% Start in on doing things  %%%%%%%%%%%%%%%%
fid = fopen(list_file, 'r');
res = textscan(fid, '%s [%d %d]');
fnames = res{1};
starts = res{2};
ends = res{3};

%select the subset of files that are desired
if isempty(file_range) file_range = 1:length(fnames); end
fnames = fnames(file_range);
starts = starts(file_range);
ends = ends(file_range);

% Video properties
%vid = struct('fname', '20121209/9085_2012-12-09-165014-0000.avi', 'timeRange', [5 250]);
%vid(2) = struct('fname', '20121209/9085_2012-12-09-171006-0000.avi', 'timeRange', [8 250]);
nfiles = length(fnames);
rew_prop = NaN*zeros(nfiles,1);
dist_prop = NaN*zeros(nfiles,1);
rew_dists = cell(nfiles,1);
distract_dists = cell(nfiles,1);
rew_dists_from_trail = cell(nfiles,1);
distract_dists_from_trail = cell(nfiles,1);
rew_dists_from_trail_persect = cell(nfiles,1);
distract_dists_from_trail_persect = cell(nfiles,1);
total_frames = NaN*zeros(nfiles,1);
frame_rate = NaN*zeros(nfiles,1);
total_turning = NaN*zeros(nfiles,1);


for ii = 1:nfiles
    fullname = fullfile(base_fname, fnames{ii});
    
    mt = class_func(fullname, [],[starts(ii) ends(ii)]);
%    mt.plotFollowing([], following_thresh, 0);
    rew_prop(ii) = mt.propTimeOnTrail([],1,following_thresh);
    dist_prop(ii) = mt.propTimeOnTrail([],2,following_thresh);
    rew_dists{ii} = mt.distanceOnTrail([],1,following_thresh);
    distract_dists{ii} = mt.distanceOnTrail([],2,following_thresh);
    rew_dists_from_trail{ii} = mt.orthogonalDistFromTrail(1:mt.nFrames,1);
    distract_dists_from_trail{ii} = mt.orthogonalDistFromTrail(1:mt.nFrames,2);
    rew_dists_from_trail_persect{ii} = mt.orthogonalDistFromTrailPerSection(1:mt.nFrames,1, following_thresh);
    distract_dists_from_trail_persect{ii} = mt.orthogonalDistFromTrailPerSection(1:mt.nFrames,2, following_thresh);
    total_turning(ii) = mt.totalTurning(1:mt.nFrames);
    % now collect a few factors about each of the movies
    total_frames(ii) = mt.nFrames;
    frame_rate(ii) = mt.frameRate;
end

% Results
ret.rew_prop = rew_prop;
ret.dist_prop = dist_prop;
ret.rew_dists = rew_dists;
ret.distract_dists = distract_dists;
ret.total_frames = total_frames;
ret.frame_rate = frame_rate;
ret.rew_dists_from_trail = rew_dists_from_trail;
ret.distract_dists_from_trail = distract_dists_from_trail;
ret.total_turning = total_turning;
ret.rew_dists_from_trail_persect = rew_dists_from_trail_persect;
ret.distract_dists_from_trail_persect = distract_dists_from_trail_persect;