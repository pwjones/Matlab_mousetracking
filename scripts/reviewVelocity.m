function vids = fixVelocities(folder_name, tracker, varargin)
% function vids = processVideoFolder(folder_name, tracker)
%
% Processes a folder of video files using the type of tracker object that is given as TRACKER in the
% input.  VARARGIN is used to specify whether or not to save the files - you can call this to load a
% group of objects in a folder easily, and not saving speeds that operation up.  Obviously, the default
% is to have the the processed tracker objects save themselves.

if ~isempty(varargin)
    saveFlag = varargin{1};
else 
    saveFlag = 1;
end

global VIDEO_ROOT;
folder_path = [VIDEO_ROOT filesep folder_name];
matfn = dir([folder_path filesep '*.mat']);
if ~isempty(matfn)
    
    s = matlabpool('size');
    if s==0 %check if parallel toolbox is running.  If not, just do regular for loop
        for ii = 1:length(matfn)
            fixVelInFile(
            basename = strtok(matfn.name(ii), '.');
            vidname = [basename '.avi']; %downside: skips other types of videos
            
            vids(ii) = trackAndSave(tracker, folder_path, fnames{ii}, starts, ends, ii, saveFlag);
        end
    else
        parfor ii = 1:length(matfn)
            vids(ii) = trackAndSave(tracker, folder_path, fnames{ii}, starts, ends, ii, saveFlag);
        end
    end
end

function mt = fixVelInFile(pathstr, filename)
