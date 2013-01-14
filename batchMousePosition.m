%let's read the file telling us what sections to load
%fid = fopen('short_video_list.txt', 'r');
fid = fopen('video_list.txt', 'r');
res = textscan(fid, '%s [%d %d]');
fnames = res{1};
starts = res{2};
ends = res{3};

base_fname = '/Users/pwjones/Movies/mouse_training/';

% Now, need to batch process these files now, detect the mouse positions
% rew_dists = cell(length(fnames)); distract_dists = cell(length(fnames));
rew_percent = zeros(length(fnames),1);  distract_percent = zeros(length(fnames),1);

for ii = 1:length(fnames)
    fullname = fullfile(base_fname, fnames{ii});
    
    mt = MouseTrackerUnder2(fullname, [],[starts(ii) ends(ii)]);
    mt.mousePosition([]);
    mt.save();
    %rew_dists{ii} = mt.distanceOnTrail([],1,5);
    rew_percent(ii) = mt.propTimeOnTrail(1:mt.nFrames,1,5);
    %distract_dists{ii} = mt.distanceOnTrail([],2,5);
    distract_percent(ii) = mt.propTimeOnTrail(1:mt.nFrames,2,5);
end