function [videoFN, folder_nums] = listBehavioralVideos(base_folder, folders, mouse_names) 
% function listBehavioralVideos(base_folder, folders, mouse_names) 
%
% Some code to read text files in each video folder and put together the list of videos for each mouse. 
% The way I see it, this way there is only need for the list of names and folders, along with the 
% text file in each individual folder, which is already used for initial processing. 
% Output: 1) a cell array of strings with the full path for each file 2) The number of folder each file 
% comes from. This is used possibly to group files by folders (days).  

%folder_nums = [];
for ii = 1:length(mouse_names)
    fni = 1; fn = {}; folderlist = [];
    for jj = 1:length(folders)
        vid_list_file = [base_folder filesep folders{jj} filesep 'tracking_times.txt'];
        if exist(vid_list_file, 'file')
            fid = fopen(vid_list_file, 'r');
        else
            fid = [];
        end
        mouse_inds = [];
        if ~isempty(fid) %check to make sure the file is open, etc
            txt = textscan(fid, '%s %d %d', 'Delimiter', ' []','MultipleDelimsAsOne', 1);
            filenames = txt{1};
            pieces = strtok(filenames, '_.');
            mouse_inds = find(strcmp(mouse_names{ii}, pieces));
        end
        for kk = 1:length(mouse_inds)
            fn{fni} = [folders{jj} filesep filenames{mouse_inds(kk)}];
            folderlist(fni) = jj;
            fni = fni+1;
            %folderlist = folders{jj};
        end
    end
    videoFN{ii} = fn;
    folder_nums{ii} = folderlist;
end
            
