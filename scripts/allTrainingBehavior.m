% allTrainingBehavior - Load up all of the mouse behavior into a big structure array 

allMiceList; % this script contains the mouse and file names to be included in the analysis.
%spring14CohortList;
base_folder = VIDEO_ROOT;
following_thresh = 20; %mm
mm_conv = .862; %mm/px
skeletonize = 0; %option to skeletonize trails or not when loading the data
clear perMouseData;
%fclose('all');

[videoList, folder_nums, day_nums] = listBehavioralVideos(base_folder, folders, mouse_names);

s = matlabpool('size');
if s~=0
    parfor ii = 1:length(videoList)
        perMouseData(ii) = loadMouseTracking(videoList{ii}, 1:length(videoList{ii}), following_thresh/mm_conv, skeletonize);
    end
else
    for ii = 1:length(videoList)
        perMouseData(ii) = loadMouseTracking(videoList{ii}, 1:length(videoList{ii}), following_thresh/mm_conv, skeletonize);
    end
end

%% Analysis segments

plotFollowingEfficiancyByDay;
xlabel('Day', 'FontSize', 18); 
ylabel('Following Efficiency (mm/sec)', 'FontSize', 18);
set(gca, 'FontSize', 18);

%% Plot the bias 

plotExpPosCDF;