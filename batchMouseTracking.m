% batch process video files for tracking
function batchMouseTracking(videoDir, movieExt, queue)

files = dir(fullfile(videoDir, ['*.' movieExt])); %selects the files with the given extension
isOpen = matlabpool('size') > 0;
if ~isOpen %start the parallel computation capabilities
    matlabpool open local_memlimited
end
% use a queue file to get limits for what sections of videos to load
if ~isempty(queue)
    queuefn = fullfile(videoDir, queue);
    qf = fopen(queuefn, 'r');
    fields = textscan(qf, '%s %f %f'); %the format for the queue file is filename_base lower_time_limit upper_time_limit
    if length(fields) == 3
        base_names = fields{1};
        lims = [fields{2} fields{3}];
    else
        disp('The queue file is not in the correct format');
    end
end

for ii=1:length(base_names)
    fn = fullfile(videoDir, [base_names{ii}, '.', movieExt]);
    if exist(fn, 'file')
        mt = MouseTracker(fn, [], lims(ii,:));
        mt.mousePosition([]);
        mt.save();
        clear mt;
    else
        disp(['The file ' fn ' does not exist. Cannot process it.']);
    end
end

matlabpool close