% Plot dataset where we've had the mice follow trails of various concentrations

concentrationDataList; % this script contains the mouse and file names to be included in the analysis.
%spring14CohortList;
base_folder = VIDEO_ROOT;
following_thresh = 20; %mm
clear perMouseData;

videoList = listBehavioralVideos(base_folder, folders, mouse_names);

s = matlabpool('size');
if s~=0
    parfor ii = 1:length(videoList)
        perMouseData(ii) = loadMouseTracking(videoList{ii}, 1:length(videoList{ii}), following_thresh);
    end
else
    for ii = 1:length(videoList)
        perMouseData(ii) = loadMouseTracking(videoList{ii}, 1:length(videoList{ii}), following_thresh);
    end
end

