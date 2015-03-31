function trackingStruct = readCCVLog(filename)
% function trackingStruct = readCCVLog(filename)
%
%
    %filename = 'trackingLog_2014-08-14-191213.txt';

    logFile = fopen(filename, 'r');
    jj = 1;
    frame = struct('areas', {}, 't', {}, 'nAreas', {});
    result = readHeader(logFile);
    while ~feof(logFile)
       segarray = textscan(logFile, 'Areas: %d');  %handle empty areas 
       for kk = 1:length(segarray{1})
           nareas = segarray{1}(kk);
           areas = struct('t', {}, 'Centroid',{}, 'Area', {}, 'PixelList', {});
           for ii=1:nareas
               tempArea = [];
               segarray = textscan(logFile, 't:%f Pos:[%f,%f] Area:%d Points:');
               tempArea.t = segarray{1};
               tempArea.Centroid = [segarray{2}, segarray{3}] + 1; % matlab pixel locations start with 1
               tempArea.Area = segarray{4};
               segarray = textscan(logFile, '[%d,%d]');
               tempArea.PixelList = [segarray{1}(:), segarray{2}(:)] + 1;
               areas(ii) = tempArea;
           end
           frame(jj).areas = areas;
           frame(jj).nAreas = length(areas);
           if (nareas > 0)
                frame(jj).t = areas(1).t;
           else
               frame(jj).t = NaN;
           end 
           jj = jj+1;
       end
    end
    fclose(logFile);
    % make sure that everything 
    trackingStruct.nPaths = result.nPaths;
    trackingStruct.paths = result.paths;
    trackingStruct.nSkelPaths = result.nSkelPaths;
    trackingStruct.skelPaths = result.skelPaths;
    trackingStruct.movieSubsample = result.movieSubsample;
    trackingStruct.binMovie = makeMovieFromAreas(frame);
    trackingStruct.frameInfo = frame;
end
% --------------------------------------------------------
function binMovie = makeMovieFromAreas(frames)
    width = 1280;
    height = 1024;
    MAX_FRAMES = 500;
    nFrames = min(MAX_FRAMES, length(frames));
    binMovie = false(height, width, nFrames);
    for i = 1:nFrames
        binMovie(:,:,i) = drawFrame(frames(i));
    end
end
% --------------------------------------------------------
function header = readHeader(logFile)
    ret = textscan(logFile, 'MovieSubsample:%d');
    header.movieSubsample = ret{1};
    ret = textscan(logFile, 'Paths:%d\n');
    header.nPaths = ret{1};
    paths = [];
    for ii = 1:header.nPaths
        ret = textscan(logFile, 'Path %u:');
        ret = textscan(logFile, '[%f,%f]');
        paths(ii).PixelList = [ret{1}, ret{2}];
        paths(ii).Area = size(paths(ii).PixelList, 1);
    end
    header.paths = paths;
    
    ret = textscan(logFile, 'skeletonPaths:%d\n');
    header.nSkelPaths = ret{1};
    for ii = 1:header.nSkelPaths
        ret = textscan(logFile, 'Path %d:');
        ret = textscan(logFile, '[%f,%f]');
        paths(ii).PixelList = [ret{1}, ret{2}];
        paths(ii).Area = size(paths(ii).PixelList, 1);
    end
    header.skelPaths = paths;
end