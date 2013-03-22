function batchMousePosition(filename, range)

fid = fopen(filename, 'r');
res = textscan(fid, '%s [%d %d]');
fnames = res{1};
starts = res{2};
ends = res{3};

base_fname = '/Users/pwjones/Movies/mouse_training/';

if isempty(range)
   range = 1:length(fnames);
elseif (range(end) > length(fnames))
   range = range(1):length(fnames);
end
fnames = fnames(range);
starts = starts(range);
ends = ends(range);

% Now, need to batch process these files now, detect the mouse positions
% rew_dists = cell(length(fnames)); distract_dists = cell(length(fnames));
rew_percent = zeros(length(fnames),1);  distract_percent = zeros(length(fnames),1);

for ii = 1:length(fnames)
    fullname = fullfile(base_fname, fnames{ii});
    
    mt = MouseTrackerUnder2(fullname, [],[starts(ii) ends(ii)]);
    mt.mousePosition([]);
    %mt.refineTracking(1:mt.nFrames);
    mt.save();
    %rew_dists{ii} = mt.distanceOnTrail([],1,5);
    rew_percent(ii) = mt.propTimeOnTrail(1:mt.nFrames,1,5);
    %distract_dists{ii} = mt.distanceOnTrail([],2,5);
    distract_percent(ii) = mt.propTimeOnTrail(1:mt.nFrames,2,5);
end