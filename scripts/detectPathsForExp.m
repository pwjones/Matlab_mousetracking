% Script to detect the paths for a folder worth of processed videos
% 10/1/2013

if ~exist('vids', 'var') %if the vids variable doesn't exist, process/load a selected folder
    base_path = VIDEO_ROOT;
    PathName = uigetdir(base_path, 'Select the folder to process');
    folders = textscan(PathName, '%s', 'Delimiter', filesep);
    folders = folders{1};
    exp_name = folders{end};
    vids = processVideoFolder(exp_name, @MouseTrackerKF);
end
vidi = 4:length(vids);
%vidi = 12:16;
%vidi = 15:length(vids);
for ii = vidi
    vids(ii).detectRefinePaths(1,1,1);
    vids(ii).save();
end