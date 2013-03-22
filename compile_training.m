function compile_training(filename, range)
%let's read the file telling us what sections to load

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

for ii = 1:length(fnames)
    fullname = fullfile(base_fname, fnames{ii});
    
    mt = MouseTrackerUnder2(fullname, [],[starts(ii) ends(ii)]);
    mt.detectPaths([],1,1);
    mt.refinePaths(1);
    mt.refinePaths(2);
    mt.save();
end

% Now, need to batch process the 