function ret = loadMouseTracking(files, file_range, following_thresh)
% function loadMouseTracking(list_file, file_range)
% 
% Function to load a set of processed videos, their saved MouseTracker
% objects. For flexibility, FILES can either be a filename or a cell array of strings 
% with individual movie fileanmes. Subsets can be selected using the FILE_RANGE argument. 
global VIDEO_ROOT;
base_fname = VIDEO_ROOT;
class_func = @MouseTrackerKF;
%class_func = @MouseTrackerUnder2;
mm_conv = 1.16; %mm/px linear
%%%%%%%%%%%%%%%% Start in on doing things  %%%%%%%%%%%%%%%%
if ~iscell(files)
    fid = fopen(files, 'r');
    res = textscan(fid, '%s [%d %d]');
    fnames = res{1};
    starts = res{2};
    ends = res{3};
    fclose(fid);
else
    fnames = files;
end

%select the subset of files that are desired
if isempty(file_range) file_range = 1:length(fnames); end
fnames = fnames(file_range);
%starts = starts(file_range);
%ends = ends(file_range);

% Video properties
%vid = struct('fname', '20121209/9085_2012-12-09-165014-0000.avi', 'timeRange', [5 250]);
%vid(2) = struct('fname', '20121209/9085_2012-12-09-171006-0000.avi', 'timeRange', [8 250]);
nfiles = length(fnames);
rew_prop = NaN*zeros(nfiles,1);
dist_prop = NaN*zeros(nfiles,1);
fake_prop = NaN*zeros(nfiles,1);
rew_dists = cell(nfiles,1);
distract_dists = cell(nfiles,1);
rew_dists_from_trail = cell(nfiles,1);
distract_dists_from_trail = cell(nfiles,1);
rew_dists_from_trail_persect = cell(nfiles,1);
distract_dists_from_trail_persect = cell(nfiles,1);
nose_trajectories = cell(nfiles,1);
traj_dir = cell(nfiles, 1);
traj_window = [];
total_frames = NaN*zeros(nfiles,1);
frame_rate = NaN*zeros(nfiles,1);
total_turning = NaN*zeros(nfiles,1);
rew_trail_area = NaN*zeros(nfiles,1);
dist_trail_area = NaN*zeros(nfiles,1);
rew_propFollowed = NaN*zeros(nfiles,1);
dist_propFollowed = NaN*zeros(nfiles,1);
dir_nums = NaN*zeros(nfiles,1);
prev_dir_name = ''; curr_dir_num = 0;

for ii = 1:nfiles
    fullname = fullfile(base_fname, fnames{ii});
    % want to number the days in case we plot performance over days
    dir_name = strtok(fnames{ii}, '/');
    if ~strcmp(dir_name, prev_dir_name);
        curr_dir_num = curr_dir_num + 1;
    end
    prev_dir_name = dir_name;
    dir_nums(ii) = curr_dir_num;
    % so this assumes that the object already exists, the video is already processed,
    % and that we are calling the constructor function to load it. 
    %mt = class_func(fullname, [],[starts(ii) ends(ii)]); 
    mt = class_func(fullname, [],[]); 
%    mt.plotFollowing([], following_thresh, 0);
    rew_prop(ii) = mt.propTimeOnTrail([],1,following_thresh);
    dist_prop(ii) = mt.propTimeOnTrail([],2,following_thresh);
    fake_prop(ii) = mt.propTimeOnTrail([], 3,following_thresh);
    rew_dists{ii} = mt.distanceOnTrail([],1,following_thresh) * mm_conv;
    distract_dists{ii} = mt.distanceOnTrail([],2,following_thresh) * mm_conv;
    rew_dists_from_trail{ii} = mt.orthogonalDistFromTrail(1:mt.nFrames,1) * mm_conv;
    distract_dists_from_trail{ii} = mt.orthogonalDistFromTrail(1:mt.nFrames,2) * mm_conv;
    rew_dists_from_trail_persect{ii} = mt.orthogonalDistFromTrailPerSection(1:mt.nFrames,1, following_thresh) * mm_conv;
    distract_dists_from_trail_persect{ii} = mt.orthogonalDistFromTrailPerSection(1:mt.nFrames,2, following_thresh) * mm_conv;
    total_turning(ii) = mt.totalTurning(1:mt.nFrames);
    [nose_trajectories{ii}, traj_dir{ii}, traj_window] = noseTrajectories(mt);
    % now collect a few factors about each of the movies
    total_frames(ii) = mt.nFrames;
    frame_rate(ii) = mt.frameRate;
    rew_trail_area(ii) = mt.paths(1).Area * mm_conv^2;
    dist_trail_area(ii) = mt.paths(2).Area * mm_conv^2;
    rew_propFollowed(ii) = mt.propTrailFollowed([], 1, 10);
    dist_propFollowed(ii) = mt.propTrailFollowed([],2,10);
end

% Results
ret.rew_prop = rew_prop;
ret.dist_prop = dist_prop;
ret.fake_prop = fake_prop;
ret.rew_dists = rew_dists;
ret.distract_dists = distract_dists;
ret.total_frames = total_frames;
ret.frame_rate = frame_rate;
ret.rew_dists_from_trail = rew_dists_from_trail;
ret.distract_dists_from_trail = distract_dists_from_trail;
ret.total_turning = total_turning;
ret.rew_dists_from_trail_persect = rew_dists_from_trail_persect;
ret.distract_dists_from_trail_persect = distract_dists_from_trail_persect;
ret.rew_trail_area = rew_trail_area;
ret.dist_trail_area = dist_trail_area;
ret.rew_propFollowed = rew_propFollowed;
ret.dist_propFollowed = dist_propFollowed;
ret.dir_nums = dir_nums;
ret.nose_trajectories = nose_trajectories;
ret.traj_window = traj_window;
ret.traj_dir = traj_dir;