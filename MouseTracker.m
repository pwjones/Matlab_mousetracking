classdef MouseTracker < handle 
    properties
        MOUSE_SCALE = 10; % minimal size of the mouse used to exclude smaller objects - radius of mask
        MIN_SIZE_THRESH = 20; %minimal size of a binary area in pixels
        MAX_SIZE_THRESH = 1200;
        videoFN = []; % filename of the video
        readerObj = ''; %videoReader object 
        framesPerSeg = 150 %the # of frames read at one time
        avgSubsample = 60 %sample every X frames to make the average
        nativeWidth = 1280 % movie dimensions, defaults
        nativeHeight = 1024
        width % after cropping
        height
        crop = [0 0 0 0]  % number of pixels trimmed from each side of movie
        frameRate %frames per second
        totalDuration %total length of the movie, sec
        nFrames
        trustFrameCount %mmread sometimes fails to exactly get the number of frames
        frameRange % frame numbers loaded in reference to original file
        avgFrame % the time averaged frame
        maxBlobs = 10 % number of mouse sections tracked
        nblobs %numer of sections found in each frame
        tailblob %the index of the tail blob for each frame
        noseblob %the index of the nose blob for each frame
        COM % center of mass for each blob
        bodyCOM %center of mass of all blobs combined
        bodyOrient %the orientation of the entire body (minus tail), in radians 0-360
        orient % orientation of each blob
        vel % velocity of each blob
        bodyVel %px/frame
        noseVel %px/frame
        direction %movement direction, blob-wise
        nosePos %the tip of the nose
        times % frame times
        areas % mouse regions
        tailVisible % boolean
        paths %binary areas for the detected paths
        boostContrast = 1; %video processing option, boolean switch to boost contrast of the movie on each frame
        % The list of properties detected in the binary image of each frame
        regprops = {'Area', 'Centroid', 'BoundingBox','MajorAxisLength','MinorAxisLength','Orientation','Extrema','PixelIdxList','PixelList','Perimeter'};
        blobsToDelete = []; %list for debugging purposes of those blobs to delete - to mark them.
        exitMovie = 0
        % Frame counter variables: There is a dedicated small area for a periodic LED signal driven by an external
        % trigger (eg physiology acquision system) in order to verify that the video is properly synced, and if it 
        % has skipped frames, to correct the error post hoc. If this is not in use, set the
        % fcArea to empty and it'll not be utilized.
        %fcArea = [1217, 959, 1262, 999]; %position in frame, [left, top, right, bottom]
        fcArea = [1206, 954, 1248, 999]
        fcLum = []; % The actual values of the LED area - sum over the area used, normalized.
        fcPeriod = 60; %the period of the repeating pattern, enables easy analysis of inter cyle variability
        syncInd = []; %this is a vector of indices for each frame corresponding to an external trigger
        % Want to assign an ID to each identified area
        blob_num = 1;
        blobID = [];
        % CCV does online tracking and writes a log of tracked areas that need to be stored and integrated into 
        % the data structures in order to analyze those paths in the same way as native tracked data. 
        logFile = 0; % Accompanying log file for the video exists
        logFN = ''; % Its name
        logStruct; % structure used to save external tracking info
    end
    properties (SetAccess = private)
        
    end
    
    methods
        function this = MouseTracker(varargin)
            % MouseTracker(varargin) - This is an object that loads a video and then detects a mouse in that video
            % letting you get the mouse position in various parts of the video.  It does this on an request basis, then 
            % saves the tracking results to minimize memory load. Also, this code is set up to track the position of
            % multiple portions of a binary thresholded image. This is meant to be different parts of the mouse when
            % viewed from below (feet, nose, and tail). Thus, the statistics compiled assume that the whole set of
            % "blobs" represents a single mouse.
            % 
            % Args: 1 - A filename or folder of the movie to load.  If there isn't one given, then a dialog is
            % presented.
            % 2 - frame range - a vector of frame numbers to consider as the whole movie, example 1:100
            % 3 - time range - a beginning and end time to consider as the whole movie, in sec, example [60 120].
            % 4 - crop
            % 5 - filename of tracked data log file that corresponds to the movie
            movie_folder = ''; filename = ''; loadVideo = 0;
            spec_frames = []; spec_timerange = []; spec_crop = [];
            if nargin > 0
                filename = varargin{1};
            end
            if nargin > 1
                spec_frames = varargin{2};
            end
            if nargin > 2
                spec_timerange = varargin{3};
            end
            if nargin > 3
                spec_crop = varargin{4};
            end
            % Figure out the filename
            if ~isempty(filename)
                if exist(filename, 'dir') % a folder was specified
                    movie_folder = filename;
                    filename = '';
                elseif exist(filename, 'file') % a valid filename is specified
                    [base_dir, base_fn, ext] =  fileparts(filename);
                    if strcmp(ext, '.mat')
                        loadVideo = 0;
                    elseif strcmp(ext, '.avi') || strcmp(ext,'.m4v') || strcmp(ext,'.mpg')
                        loadVideo = 1;
                    elseif strcmp(ext, '.txt')
                        this.logFile = 1;
                        this.logFN = filename;
                        loadVideo = 0;
                    end
                else
                    loadVideo = 0;
                end
            else
                loadVideo = 1;
            end
            % Trying to figure out how to load, and whether prompt for a filename
            if (isempty(filename) || ~exist(filename)) && loadVideo % && ~mclIsNoDisplaySet() %a non-complete path was specified, and needs to be chosen
                movie_folder = '/Volumes/Alexandria/pwj_data/mouse_training/';
                [filename, movie_folder] = uigetfile([movie_folder '*.*']);
                [~, base_fn, ext] =  fileparts(filename);
                base_dir = movie_folder;
                if filename == 0 % the user has canceled the file selection
                    return;
                end
                this.videoFN = [movie_folder filename];
            elseif ~loadVideo
                this.videoFN = '';
            else
                this.videoFN = filename;
            end
            
            if ~isempty(spec_crop) %there is a crop vector present
                this.crop = varargin{4};
                if isempty(this.crop)
                    this.crop = [0 0 0 0];
                end
                this.width = this.readerObj.width - sum(this.crop(1:2));
                this.height = this.readerObj.height - sum(this.crop(3:4));
            else
                this.crop = [ 0 0 0 0];
                this.width = this.nativeWidth;
                this.height = this.nativeHeight;
            end
            
            if (nargin > 4) % there is an additional log file given
                this.logFile = 1;
                this.logFN = varargin{5};
            end
            
            mat_fn = fullfile(base_dir, [base_fn '.mat']); % make a '.mat' filename
            % Branch on what to load 1) mat file 2) video file 3) log file
            if exist(mat_fn, 'file')  % then load the saved mat file
                disp(['Loading the saved object for: ' mat_fn]);
                arg_vidFN = this.videoFN;
                try
                    saved = load(mat_fn);
                catch ME
                    disp(['Error loading: ' mat_fn]);
                end
                if ~strcmp(class(saved.this), class(this))
                    saved = eval([class(this) '(saved.this);']);
                    this = saved;
                else
                    this = saved.this;
                end
                if ~isempty(this.videoFN)
                    if exist(this.videoFN, 'file')
                        this.findMovieFile(this.videoFN);
                    elseif (exist(arg_vidFN)) %if files are moved, resolve
                        this.findMovieFile(arg_vidFN);
                    else
                        error('There is no valid movie file specified. Cannot load object');
                    end
                end
                this.videoFN = arg_vidFN; %for some reason, in no case does the full video filename gets saved
            elseif loadVideo %initialize normally, using movie 
                this = this.initWithMovie(this.videoFN, spec_timerange, spec_frames);
                if this.logFile
                    this.initWithLog(this.logFN);
                end
                %initialize tracking data
                this.clearCalcData();
            else
                this.nFrames = 0;
                if this.logFile
                    % initialize based on the log file instead of a movie
                    this = this.initWithLog(this.logFN);
                end
                this.clearCalcData();
            end
        end %function MouseTracker
        
         % ------------------------------------------------------------------------------------------------------
        function this = initWithMovie(this, filename, spec_timerange, spec_frames) 
            
            % Read just a couple of frames to get an idea of video speeds, etc.
            this.readerObj = VideoReader(this.videoFN);
            this.frameRate = this.readerObj.FrameRate;
            this.nativeWidth = this.readerObj.Width;
            this.nativeHeight = this.readerObj.Height;
            this.fcPeriod = round(this.frameRate);
            
            % figure out the range of frames to be considered
            this.frameRange = [];
            if ~isempty(spec_timerange) % figure out the frames from times
                time_range = spec_timerange;
                t_offset = time_range(1);
                fr = round((time_range)*this.frameRate + 1);
                if (fr(2) > this.readerObj.NumberOfFrames)
                    fr(2) = this.readerObj.NumberOfFrames;
                end
                if isempty(spec_frames) % use frames if specified
                    this.frameRange = fr(1):fr(2);
                else
                    this.frameRange = spec_frames;
                end
                
            end
            if isempty(this.frameRange) % then specify the whole video's range
                this.frameRange = 1:this.readerObj.NumberOfFrames; %undershoot to avoid potential problems
                this.trustFrameCount = 1;
            end
            
            % read in enough frames to create the average frame
            disp('MouseTracker: reading movie to make average frame');
            avgRange = 1:this.avgSubsample:length(this.frameRange);
            while (length(avgRange) < 20)
                this.avgSubsample = floor(this.avgSubsample/2);
                avgRange = 1:this.avgSubsample:length(this.frameRange);
            end
            vidArray = this.readFrames(this.listToIntervals(avgRange), 'discontinuous');
            this.avgFrame = uint8(mean(vidArray, 3));
            
            % Populate the fields with empty data
            this.nFrames = length(this.frameRange);
            this.times = ((1:this.nFrames)-1)/this.frameRate; % time vector starts at 0
        end
        
        % ------------------------------------------------------------------------------------------------------
        function this = initWithLog(this, filename) 
            if isempty(this.videoFN)
                this.readerObj = '';
            end
            this.readTrackingLog(filename);
            if ~isempty(this.videoFN)
                this.alignTrackingLog();
            end
            this.nFrames = length(this.logStruct.frameInfo);
            this.frameRange = 1:this.nFrames;
            this.avgFrame = zeros(this.height, this.width);
            
            % Need to construct a time vector with uneven sampling, and
            % filling in those entries when there were no detected objects
            this.times = zeros(this.nFrames, 1);
            t = [this.logStruct.frameInfo.t];
            first_t = find(t>0, 1,'first');
            dt = diff(t); 
            mean_dt = nanmean(dt);
            toff = t(first_t) - (first_t*mean_dt);
            this.times(1) = 0;
            for ii = 2:this.nFrames
                this.times(ii) = this.logStruct.frameInfo(ii).t - toff;
                if isnan(this.times(ii))
                    this.times(ii) = this.times(ii-1) + mean_dt;
                end
            end
            this.frameRate = 1000/mean_dt;
        end
        
        % ------------------------------------------------------------------------------------------------------
        function [wantedCOM, nose] = mousePosition(this, time_range)
        % function [wantedCOM, nose] = mousePosition(this, time_range)    
         % Returns the position of the mouse center of mass (body) and 
         % the nose for the specified time range
            wantedFrames = this.timesToFrames(time_range); 
            wantedCOM = this.bodyCOM(wantedFrames,:);
            if isnan(wantedCOM)
                nanFrames = wantedFrames(isnan(wantedCOM(:,1)));
                this = this.findMouse(nanFrames);
            end 
            % get what they wanted
            wantedCOM = this.bodyCOM(wantedFrames,:);
            nose = this.nosePos(wantedFrames, :);
        end
        
        % ------------------------------------------------------------------------------------------------------
        function vel = mouseVelocity(this, time_range)
            % function vel = mouseVelocity(this, time_range)
            % Returns the velocity of the mouse body for the specified frames
            frames = this.timesToFrames(time_range);
            vel = this.vel(frames);
            if isnan(vel)
                this.mousePosition(time_range);
            end
            vel = this.computeVelocity(frames);
        end
        
        % ----------------------------------------------------------------------------------------------------
        function tailPos = tailPosition(this, frames)
        % function tailPos = tailPosition(this, frames)
        %
        % We aren't going to want the tail position all that much, so assemble it each time
            tv = logical(this.tailVisible(frames));
            tb = this.tailblob(frames);
            tailPos = NaN*zeros(length(frames), 2);
            areas = this.areas(frames,:);
            for ii = 1:length(areas)
                if tv(ii)
                    tailPos(ii,:) = areas(ii,tb(ii)).Centroid;
                end
            end
        end
        
        % ------------------------------------------------------------------------------------------------------
        function ah = plotVelocity(this, frames)
            % function plotVelocity(this) 
            if ~exist('frames', 'var') || isempty(frames); frames = 1:this.nFrames; end
            figure;
            sorted_vel = sort(this.vel); %sorted vector of speeds for colormapping
%             pathIm = this.plotPaths()*255;
%             redIm = min(this.avgFrame + pathIm, 255);
%             blueIm = max(this.avgFrame - pathIm, 0);
%             bgIm = cat(3, blueIm, blueIm, redIm);
%             imshow(bgIm); hold on;
            imshow(this.plotPathsOnBg()); hold on;
            xlim([0 this.width]); %fit the axes to the image
            ylim([0 this.height]);
            %cm = colormap(spring(this.nFrames));
            vel_vect = this.vel(frames);
            [cm cinds] = getIndexedColors('jet', vel_vect, 1);
            %cm = brighten(cm, .5);
            % unfortunately, I haven't been able to figure out an easier way to do colormapping of points plotted over
            % a B&W image, without the points and the image sharing the same colormap.  So, here I'm just plotting every
            % point separately, with a color indicative of the velocity of the animal.
            for ii=1:length(frames)
                fi = frames(ii);
                %ci = find(sorted_vel == this.vel(fi), 1, 'first');
                if ~isnan(vel_vect(ii))
                    line('Xdata', this.nosePos(fi,1), 'Ydata', this.nosePos(fi,2), 'Marker', '.', 'MarkerSize', 10, 'Color', cm(cinds(ii),:));
                end
            end 
            ah = gca;
            %Make another figure entirely to get a color scale
            fh2 = figure;
            %ah2 = axes('position', [.1 .1 .7 .8], 'visible', 'off');
            colormap(cm);
            %axes(ah2);
            pcolor([this.vel(frames,1), this.vel(frames,1)]);
            %line('parent', ah2, 'ydata', sorted_dir, 'xdata', 1:length(sorted_dir)); 
            colorbar;
        end
        % ------------------------------------------------------------------------------------------------------
        function plotVelocityTimes(this, time_range)
         % Version using a range of times, rather than frames
            frames = this.timesToFrames(time_range);
            this.plotVelocity(frames);
        end
        
        % ------------------------------------------------------------------------------------------------------
        function plotPosition(this, frames, varargin)
            % function plotPosition(this, frames)
            %
            % This plots the detected mouse positions in red over the background frame
            c = [.3, .6, 1]; %plotblue
            marker = '.';
            ah = [];
            overwrite = 1;
            if length(varargin) >=1
                ah = varargin{1};
                newFig = 0;
            end
            if length(varargin) >= 2
                overwrite = varargin{2};
            end
            if length(varargin) >=3
                c = varargin{3};
            end
            if length(varargin) >=4
                marker = varargin{4};
            end
            if isempty(frames); frames = 1:this.nFrames; end
            if isempty(ah)
                fh = figure; 
                ah = axes('position', [.1, .1, .7 .8]); hold on;
                newFig = 1;
            end
            if (newFig || overwrite) 
                pathIm = this.plotPathsOnBg();
                imshow(pathIm); hold on;
            end
            xlim([0 this.width]);
            ylim([0 this.height]);
            for ii=1:length(frames)
                fi = frames(ii);
                if ~isnan(this.noseblob(fi))
                    line('Parent', ah, 'Xdata', this.nosePos(fi,1,1), 'Ydata', this.nosePos(fi,2,1), ...
                        'Marker', '.', 'MarkerSize', 6, 'Color', c);
                end
            end 
        end
        
        % ------------------------------------------------------------------------------------------------------
        function plotNosePosition(this, frames)
            % function plotDirection(this, frames)
            plotblue = [.3, .6, 1]; % a nice blue color
            if ~exist('frames', 'var') || isempty(frames); frames = 1:this.nFrames; end
            fh = figure;
            ah = axes('position', [.1, .1, .7 .8]);
            bgIm = this.plotPathsOnBg();
            imshow(bgIm); hold on;
            xlim([0 this.width]);
            ylim([0 this.height]);
            for ii=1:length(frames)
                fi = frames(ii);
                if ~isnan(this.noseblob(fi))
                    line('Parent', ah, 'Xdata', this.areas(fi,this.noseblob(fi)).Centroid(1), 'Ydata', ...
                        this.areas(fi,this.noseblob(fi)).Centroid(2), 'Marker', '.', 'MarkerSize', 8, 'Color', plotblue);
                end
            end 
            title(this.videoFN);
        end
        % ------------------------------------------------------------------------------------------    
        function [nb_stat, nose_bright] = getNoseBrightness(this, frames, operation)
            % operation is a function pointer to operation to be performed
            % on the set of nose pixels 
            nSegs = ceil(length(frames)/this.framesPerSeg);
            movieDone = 0; 
            dbg = 0;
            nose_bright = NaN*zeros(length(frames), this.MAX_SIZE_THRESH);
            nb_stat = NaN*zeros(length(frames),1);
            for jj = 1:nSegs %divide the computation up into segments so that there is never too much in memory at once
                first = (jj-1)*this.framesPerSeg + 1; %first frame of segment
                if (jj == nSegs) % the last frame
                    last = length(frames);
                else
                    last = (jj)*this.framesPerSeg;
                end
                disp(sprintf('Finding brightnesses for frames %d - %d', first, last));
                segFrames = frames(first:last);
                
                frameArray = this.readMovieSection(segFrames,'diff');
                for ii = 1:(last-first)
                    fi = frames(first+ii-1);
                    if ~isnan(this.noseblob(fi))
                        ni = this.areas(fi, this.noseblob(fi)).PixelIdxList;
                        frame = frameArray(:,:,ii);
                        nose_px = frame(ni);
                        nose_bright(first+ii-1, 1:length(nose_px)) = nose_px; 
                        nb_stat(first+ii-1) = operation(nose_px(:));
                    end
                end
                        
            end
        end
        % -------------------------------------------------------------------------------------------
        function plotDirection(this, frames)
            % function plotDirection(this, frames)
            plotblue = [.3, .6, 1]; % a nice blue color
            if ~exist('frames', 'var') || isempty(frames); frames = 1:this.nFrames; end
            fh = figure;
            ah = axes('position', [.1, .1, .7 .8]);
            cm = colormap(hsv(this.nFrames));
            sorted_dir = sort(this.direction); %sorted vector for colormapping
            bgIm = this.plotPathsOnBg();
            imshow(bgIm); hold on;
            xlim([0 this.width]);
            ylim([0 this.height]);
            % unfortunately, I haven't been able to figure out an easier way to do colormapping of points plotted over
            % a B&W image, without the points and the image sharing the same colormap.  So, here I'm just plotting every
            % point separately, with a color indicative of the velocity of the animal.
            for ii=1:length(frames)
                fi = frames(ii);
                ci = find(sorted_dir == this.direction(fi), 1, 'first');
                if ~isnan(this.noseblob(fi))
                    line('Xdata', this.areas(fi,this.noseblob(fi)).Centroid(1), 'Ydata', ...
                        this.areas(fi,this.noseblob(fi)).Centroid(2), 'Marker', '.', 'MarkerSize', 8, 'Color', cm(ci,:));
                end
                if ~isnan(ci)
                    %line('Xdata', this.COM(fi,1,1), 'Ydata', this.COM(fi,2,1), 'Marker', '.', 'MarkerSize', 10, 'Color', cm(ci,:));
                    %line('Xdata', this.COM(fi,1,1), 'Ydata', this.COM(fi,2,1), 'Marker', '.', 'MarkerSize', 10, 'Color', 'r');
                    %line('parent', ah2, 'Xdata', this.COM(ii,1), 'Ydata', this.COM(ii,2), 'Marker', '.', 'MarkerSize', 10, 'Color', cm(ci,:));
                end
            end 
            title(this.videoFN);
            % Make another graph with the color scale
            fh2 = figure;
            ah2 = axes('position', [.1 .1 .7 .8], 'visible', 'off');
            cm = colormap(hsv(this.nFrames));
            axes(ah2);
            pcolor(repmat(sorted_dir(:), [1 3]));
%             line('parent', ah2, 'ydata', sorted_dir, 'xdata', 1:length(sorted_dir)); 
            colorbar;
        end
        
        % ---------------------------------------------------------------------------
        function plotFollowing(this, frames, dist_thresh, textflag)
            % function plotFollowing(this, frames, dist_thresh, textflag)
            % plot the following segments based on the following information 
            if isempty(frames)
                frames = 1:this.nFrames;
            end
            [dists, fframes] = this.distanceOnTrail(frames,1,dist_thresh);
            this.plotNosePosition(frames); hold on;
            for i=1:size(fframes,1)
                range = fframes(i,1):fframes(i,2);
                np = this.nosePos(range, :);
                plot(np(:,1), np(:,2), '.m', 'MarkerSize',10);
                if textflag
                    text(mean(np(:,1)), mean(np(:,2)), num2str(dists(i)), 'Color', 'w', 'FontSize', 12);
                end
            end
        end
        
        % ---------------------------------------------------------------------------
        function ah = plotFollowingSide(this, frames, dist_thresh, textflag)
            % function ah = plotFollowingSide(this, frames, dist_thresh, textflag)
            % plot the following segments based on the following information
            if isempty(frames)
                frames = 1:this.nFrames;
            end
            [dists, fframes] = this.distanceOnTrail(frames,1,dist_thresh);
            this.plotNosePosition; hold on;
            for i=1:size(fframes,1)
                range = fframes(i,1):fframes(i,2);
                np = this.nosePos(range, :);
                %incorporate the side of the trail coloring
                signed_dists = this.orthogonalDistFromTrail(range,1);
                negd = signed_dists <= -1;
                posd = signed_dists >= 1;
                plot(np(posd,1), np(posd,2), '.m', 'MarkerSize',10);
                plot(np(negd,1), np(negd,2), '.y', 'MarkerSize',10);
                if textflag
                    text(mean(np(:,1)), mean(np(:,2)), num2str(dists(i)), 'Color', 'w', 'FontSize', 12);
                end
            end
            ah = gca;
        end
        
        % ---------------------------------------------------------------------------
        function plotOrientation(this, frames)
            % function plotOrientation(this, frames)
            if ~exist('frames', 'var') || isempty(frames); frames = 1:this.nFrames; end
            figure;
            cm = colormap(hsv(this.nFrames));
            sorted_orient = sort(this.bodyOrient); %sorted vector of speeds for colormapping
            imshow(this.avgFrame); hold on;
            xlim([0 this.width]);
            ylim([0 this.height]);
            % unfortunately, I haven't been able to figure out an easier way to do colormapping of points plotted over
            % a B&W image, without the points and the image sharing the same colormap.  So, here I'm just plotting every
            % point separately, with a color indicative of the velocity of the animal.
            for ii=1:length(frames)
                fi = frames(ii);
                ci = find(sorted_orient == this.bodyOrient(fi), 1, 'first');
                if ~isnan(ci)
                    %line('Xdata', this.COM(ii,1), 'Ydata', this.COM(ii,2), 'Marker', '.', 'MarkerSize', 10, 'Color', cm(ci,:));
                    line('Xdata', this.nosePos(fi,1), 'Ydata', this.nosePos(fi,2), 'Marker', '.', 'MarkerSize', 10, 'Color', cm(ci,:));
                end
            end 
            % Make another graph with the color scale
            fh2 = figure;
            ah2 = axes('position', [.1 .1 .7 .8], 'visible', 'off');
            cm = colormap(hsv(this.nFrames));
            axes(ah2);
            pcolor(repmat(sorted_orient(:), [1 3]));
            %line('parent', ah2, 'ydata', sorted_dir, 'xdata', 1:length(sorted_dir)); 
            colorbar;
        end
        
        % ------------------------------------------------------------------------------------------------------
        function writeMovie(this, filename, movieType, frames, dispCrop)
            % function writeMovie(this, filename, movieType, frames, dispCrop)
            % Writes a movie of the tracked mouse to disk at the given location
            if isempty(frames)
                frames = 1:this.nFrames;
            end
            vidWriter = VideoWriter(filename, 'Motion JPEG AVI');
            open(vidWriter);
            figure;
            for ii=frames
                this.showFrame(ii, movieType, dispCrop);
                currFrame = getframe; % now to the business of writing the movie
                writeVideo(vidWriter,currFrame);
                hold off;
            end
            close(vidWriter);
        end
        
        % -------------------------------------------------------------------------------------------------
        function showMovie(this, movieType, frames, varargin)
            % function showMovie(this, useBinMovie, frames, dispCrop[OPTIONAL], overlay[OPTIONAL])
            if nargin >=4
                dispCrop = varargin{1};
            else
                dispCrop = [];
            end
            if nargin >=5
                overlayMov = varargin{2};
            end
            fh = figure;
            this.exitMovie = 0;
            set(fh, 'WindowKeyPressFcn', @this.exitMovieLoop);
            if ~exist('frames','var') || isempty(frames)
                frames = 1:this.nFrames;
            end
            for ii=1:length(frames)
                fi = frames(ii);
                if this.exitMovie
                    break
                end
                if nargin < 5
                    this.showFrame(fi, movieType, dispCrop);
                else
                    this.showFrame(fi, movieType, dispCrop, overlayMov(:,:,ii));
                end
                hold off;
                pause(1/this.frameRate/1); %faster than the natural framerate due to impatience.
            end
            this.exitMovie = 0;
            set(fh, 'WindowKeyPressFcn', '');
        end
        
        % ------------------------------------------------------------------------------------------------------
        function exitMovieLoop(this, src, event)
            % Function to set a flag internally to exit a movie that is being displayed.
            % It is not used for any other purpose
            if strcmp(event.Key, 'q') || strcmp(event.Key, 'escape')
                this.exitMovie = 1;
            end
        end
        
        function save(this, varargin)
            % Call this function with the filename that you want to save as, otherwise the default is the
            % name of the videofile.  It will save in the video file directory.
            if ~isempty(this.videoFN)
                [basedir, fn_base, fn_ext] = fileparts(this.videoFN);
            else
                [basedir, fn_base, fn_ext] = fileparts(this.logFN);
            end
            if ~isempty(varargin)
                % save with the given filename
                fn_base = varargin{1};
            end
            full_fn = fullfile(basedir, [fn_base '.mat']);
            save(full_fn, 'this');
            disp(['Saved object in ' full_fn]);
        end
        
        function mov = returnFrames(this, frames, movieType)
            % mov = returnFrames(this, frames, binary)
            % returns the frames of the movie specified, binary or greyscale
            mov = this.readMovieSection(frames, movieType);
            
        end
        
        %------------------------------- Frame Counter functions --------------------------------------
        function setFrameCountArea(this)
        % Brings up a selector to override the default area
            figure;
            imshow(this.avgFrame)
            H = imrect;
            pos = round(H.getPosition);
            this.fcArea = [pos(1) pos(2) pos(1)+pos(3) pos(2)+pos(4)];
        end
        
        function missing = isMissingFrame(this)
        % Check if frame(s) is(are) missing
            missing = 0;
            if ~isempty(this.fcArea)
                missing = isMissingFrame(this.fcLum, this.fcPeriod, .02); %functionality outsourced
            end 
        end
        
        function updateFCLum(this)
           cycles = ceil(this.nFrames/this.framesPerSeg);
           endFrame = 0;
           for ii = 1:cycles
               disp(sprintf('Starting segment %d of %d', ii, cycles));
               beginFrame = endFrame + 1;
               endFrame = min(beginFrame + this.framesPerSeg - 1, this.nFrames);
               rawArray = this.readFrames([beginFrame endFrame], 'continuous');
               if ~isempty(this.fcArea) % excludes the area used for counting frames
                    fca = this.fcArea;
                    fc_mov = rawArray(fca(2):fca(4), fca(1):fca(3), :);
                    fc_val = squeeze(sum(sum(fc_mov,1),2));
               else
                    fc_val = NaN*ones(length(beginFrame:endFrame), 1);
               end
               this.fcLum(beginFrame:endFrame) = fc_val;
           end
        end
    %end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRIVATE METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %methods %(Access = private)
        
        function this = findMouse(this, frames)
            % function this = findMouse(this, frames)
            
            nSegs = ceil(length(frames)/this.framesPerSeg);
            movieDone = 0; 
            dbg = 0;
            
            for jj = 1:nSegs %divide the computation up into segments so that there is never too much in memory at once
                first = (jj-1)*this.framesPerSeg + 1; %first frame of segment
                if (jj == nSegs) % the last frame
                    last = length(frames);
                else
                    last = (jj)*this.framesPerSeg;
                end
                segFrames = frames(first:last);
                disp(['Finding mouse in segment ' num2str(jj)]);
                if ~this.logFile
                    frameArray = this.readMovieSection(segFrames,'bin');
                else
                    frameArray = this.makeMovieFromLog(segFrames);
                end
                segFrames = segFrames(1:size(frameArray, 3));
                newFrameCount(jj) = this.nFrames;
                trimmed(jj) = 0; 
                if size(frameArray,3) < length(segFrames) %there was a problem reading some frames, so we need to adjust
                    disp('In MouseTracker.findMouse, frameArray has too few frames, revising expectations');
                    fdiff = this.nFrames - segFrames(end); %now, consider the entire length of the frames
                    newFrameCount(jj) = this.nFrames - fdiff; 
                    trimmed(jj) = 1;
                    %movieDone = 1; %don't read movie anymore
                end
                for ii = 1:length(segFrames)  %have to work 1 frame at a time, unfortunately
                    lm = bwlabel(frameArray(:,:,ii)); %label matrix
                    fi = segFrames(ii);
                    temp_reg = regionprops(lm, this.regprops); %gives the labels for the areas
                    if ~isempty(temp_reg)
                        if (1 == 0) %debug plotting
                           bin_im = zeros(this.height, this.width); 
                           %bin_im(temp_area.PixelIdxList) = 255;
                           figure; imshow(lm, [0 max(lm(:))]); colormap(jet);
                        end
                        [~, areai] = sort([temp_reg.Area], 'descend');
                        temp_reg = temp_reg(areai); %reorder in terms of area 
                        area_sel = ([temp_reg.Area] >= this.MIN_SIZE_THRESH & [temp_reg.Area] <= this.MAX_SIZE_THRESH);
                        temp_reg= temp_reg(area_sel);
                        for kk=1:length(temp_reg)
                            temp_reg(kk).Orientation = -temp_reg(kk).Orientation; %the returned orientation is not correct, need to correct
                        end
                        this.nblobs(fi) = min(this.maxBlobs, length(temp_reg)); 
                        this.areas(segFrames(ii),1:this.nblobs(fi)) = temp_reg(1:this.nblobs(fi)); %only keep the top sized blobs
                       
                        this.detectTail(fi);
                        if dbg %debug plotting
                           bin_im = zeros(this.height, this.width); 
                           for kk=1:length(temp_reg)
                               bin_im(temp_reg(kk).PixelIdxList) = kk;
                           end
                           figure; imshow(label2rgb(bin_im,'jet','k'));
                           hold on; plot(this.COM(fi,1,this.tailblob(fi)), this.COM(fi,2,this.tailblob(fi)),'r+');
                           hold on; plot(this.bodyCOM(fi,1), this.bodyCOM(fi,2),'yx','MarkerSize',12);
                        end
                    else %this is a kluge - unclear what the best thing to do is if you don't detect a blob
                        subI = max(1, ii-1);
                        this.areas(segFrames(ii)) = this.areas(frames(subI)); %this will error on ii=1
                    end
                end
                clear frameArray;
            end
            if trimmed %means the last set of frames came out shorter than expected from movie info, trimming things
                this.nFrames = min(newFrameCount);
                frames = frames(1:this.nFrames); %just trim them off the end - this error only happens for last seg
                this.trimFields(); %trim the object property arrays
            end
            this.refineTracking(frames);    
        end
        
        % ------------------------------------------------------------------------------------------------------
        function refineTracking(this, frames)
        % function refineTracking(this, frames)
        %
        % It makes sense to group these because they should all be performed together, and recomputed if the obh.
            this.removeStaticObjects();
            this.detectTail(frames);
            this.bodyCOM(frames,:) = this.computeBodyPos(frames,0);
            this.bodyOrient(frames) = this.computeBodyOrientation(frames);
            this.bodyOrient(frames) = this.fixOrientation(frames); %this should straighten out the orientations
            this.nosePos(frames, :) = this.findNose(frames); 
            this.vel(frames, :) = this.computeVelocity(frames);
        end
        
        % -----------------------------------------------------------------------------------------------------
        function detectTail(this, frames)
            % Trying to detect if we can see a tail on a frame by frame basis
            for ii=1:length(frames)
                regions = this.areas(frames(ii),:);
                % going to find the blob with the highest eccentricity to be the tail
                ecc = [regions.MajorAxisLength]./[regions.MinorAxisLength];
                [max_ecc, maxi] = nanmax(ecc);
                this.tailblob(frames(ii)) = maxi;
                if max_ecc>3
                    this.tailVisible(frames(ii)) = 1;
                else
                    this.tailVisible(frames(ii)) = 0;
                end
            end
        end
        % ------------------------------------------------------------------------------------------------------
        function nosep = findNose(this, frames)
        % function nosep = findNose(this, frames, blobi)
        %
        % So, the method for finding which blob is the nose is to take the tail and body center and take the 
        % distance along that vector for every blob centroid. The one with the largest is the nose.
            nosep = NaN* zeros(length(frames), 2);
            for ii = 1:length(frames)
                fi = frames(ii);
                if (this.nblobs(fi) > 3) %only works if there are enough things to track
                    if this.tailVisible(fi) && (this.nblobs(fi) > 2)%we can only do this first method if the tail is present
                        % compute the vector from tail to body center
                        bodyVect = this.bodyCOM(fi,:) - this.areas(fi,this.tailblob(fi)).Centroid;
                        dist = [];
                        for jj=1:this.nblobs(fi)
                            areaVect = this.areas(fi,jj).Centroid - this.areas(fi,this.tailblob(fi)).Centroid;
                            dist(jj) = dot(areaVect, bodyVect);
                        end
                        [maxd, maxi] = max(dist); %the maximum distance and ind along tail-body vector of each blob  
                        this.noseblob(fi) = maxi; 
                        nosep(ii,:) = this.areas(fi,this.noseblob(fi)).Centroid; % The nose is the farthest along vector
                    end
                end
            end
            % putting in a maximum distance criterion of 500 px away from last frame 
            dists = diff(nosep,1,1);
            dists = sqrt(sum(dists.^2,2));
            toofar = dists > 100;
            nosep(toofar,:) = NaN*zeros(sum(toofar), 2);
        end
        
        % ------------------------------------------------------------------------------------------------------
        function assignBlobIDs(this, frame)
            % Assigns unique IDs to blobs based on their overlap with past blobs.  If you have a >30% overlap, you are assumed to 
            % be the same blob and therefore get teh same identifier.  If you don't, you get a new one.  
            % Some parameters for assignment:
            FTF_OL = .2; % The frame-to-frame overlap to be called the same area
            LOW_OL = .1; % A lower overlap threshold for determining merging
            MAX_DIST = ceil(500/1.16/this.frameRate);% The maximum distance to spots can be to autmatically be given same ID.
            MAX_SIZE_DIFF = .3;
            
            fi = frame;
            if (frame ~= 1)
                prev = this.areas(frame-1, 1:this.nblobs(fi-1));
                matched = zeros(this.nblobs(frame),1); %does each blob have a prev frame match - if so which
                for ii = 1:this.nblobs(fi)
                    blob = this.areas(frame, ii);
                    if ~isempty(prev) %if there were blobs on the last frame, checked to avoid error
                        overlap = regionOverlap(blob, prev, this.regprops);
                        oli = find(overlap >= LOW_OL);
                        ol_vals = overlap(oli); %overlap values
                        if length(oli > 1) % Sort in descending overlapping-ness
                            [ol_vals, sorti] = sort(ol_vals(:), 1, 'descend');
                            oli = oli(sorti);
                        end
                        prev_pos = reshape([prev(oli).Centroid], 2,[])';
                        com_dist = NaN*zeros(size(prev_pos,1),2);
                        for jj=1:size(prev_pos,1) 
                            com_dist(jj,:) = blob.Centroid - prev_pos(jj,:); 
                        end
                        com_dist = sqrt(sum(com_dist.^2, 2));
                        npx = [prev(oli).Area];
                        size_diff = abs(npx - blob.Area)/blob.Area;
                        % loop through overlapping blobs in prev frames, looking for the best qualified one
                        % that isn't taken.
                        for jj = 1:length(oli)
                            %test if they are similar enough to be considered same
                            size_good = (size_diff(jj) < MAX_SIZE_DIFF);
                            dist_good = com_dist(jj) < MAX_DIST;
                            overlap_good = ol_vals(jj) > FTF_OL;
                            if ((dist_good + overlap_good) >= 1) && size_good
                                if sum(matched == oli(jj))<1 % is that ID already in this frame?
                                    matched(ii) = oli(jj);
                                end
                            end
                        end
                    end
                    if matched(ii)
                        this.blobID(frame, ii) = this.blobID(frame-1, matched(ii));
                    else 
                        this.blobID(frame, ii) = this.blob_num;
                        this.blob_num = this.blob_num+1;
                    end
                end
            else
                for ii = 1:this.nblobs(fi)
                    this.blobID(frame, ii) = this.blob_num;
                    this.blob_num = this.blob_num+1;
                end
            end
        end
        
        % ------------------------------------------------------------------------------------------------------
        function bodyCOM = computeBodyPos(this, frames, includeTail)
            %
            if isempty(includeTail) includeTail=0; end
            bodyCOM = NaN*zeros(length(frames), 2);
            for i = 1:length(frames)
                fi = frames(i);
                taili = this.tailblob(fi);
                if (isnan(taili)) taili = this.maxBlobs + 1; end
                temp_pos = [];
                temp_areas = this.areas(fi,1:this.nblobs(fi));
                for j=1:this.nblobs(fi)
                    this.orient(fi, j) = [this.areas(fi,j).Orientation]./ 180 * pi; %this is the rough estimate
                    positions = permute(reshape([temp_areas(j).Centroid], 2,[]), [3 1 2]);
                    this.COM(fi, :, j) = positions;
                    if (includeTail) || (j ~= taili) %exclude the tail in the body position calculation
                        temp_pos = cat(1, temp_pos, temp_areas(j).PixelList); %for making the mean weighted by size
                        %temp_pos = cat(1, temp_pos, positions); %making it not weighted by size
                    end
                end
                bodyCOM(i, :) = mean(temp_pos);
            end
        end
        
       
        % ------------------------------------------------------------------------------------------------------
        function orient = computeBodyOrientation(this, frames)
        % function [this, orient] = computeBodyOrientation(this, frames)
        %
        % The idea here is to fill in the body space using dilation of the individual body blobs (excluding tail)
        % in order to get a reading of the animal orientation.
            se = strel('disk',20);
            orient = NaN*zeros(length(frames),1);
            for ii=1:length(frames)
                fi = frames(ii);
                if (this.nblobs(fi) > 1)
                    areas = this.areas(fi,1:this.nblobs(fi));
                    if (this.tailVisible(fi))
                        sel = true(this.nblobs(fi), 1);
                        sel(this.tailblob(fi)) = 0;
                        areas = areas(sel);
                    end
                    bw = false(this.height, this.width);
                    for jj=1:length(areas)
                        bw(areas(jj).PixelIdxList) = 1;
                    end
                    bw2 = imdilate(bw,se);
                    propstr = {'Orientation'};
                    props = regionprops(bw2, propstr);
                    if (~isempty(props))
                        orient(ii) = props(1).Orientation ./ 180 .* pi;
                    else 
                        orient(ii) = NaN;
                    end
                end
            end
            %this.orient(frames) = orient;
        end
        
        % ------------------------------------------------------------------------------------------------------
        function vel = computeVelocity(this, frames)
            % function vel = getVelocity(this, frames)
            % Computes the velocity as the frame by frame difference in position
            
            %in order to get a velocity for the first frame (0,0) position is assumed at frame 0 
            this.vel = NaN*zeros(size(this.bodyCOM,1),size(this.bodyCOM,3));
            vel = NaN*zeros(length(frames), size(this.bodyCOM,3));
            COM = this.computeBodyPos(frames, 0); %try including the tail in the body position for velocity calculation
            %COM = this.nosePos(frames,:);
            diff_com = diff(COM, 1, 1); %differences in COM for each x,y position
            dists = [sqrt(sum(COM(1,:).^2, 2)); sqrt(sum(diff_com.^2, 2))]; %euclidian distances of each frame transition
            %this.vel(frames) = dists * this.frameRate; % Velocities in px/sec
            this.vel(frames) = dists;
            vel = this.vel(frames);
            diff_com = [COM(1,:); diff_com]; % add the first position to get an equal sized vector
            this.direction(frames) = cart2pol(diff_com(:,1), diff_com(:,2)); %this is the direction of motion
        end
        
        % ------------------------------------------------------------------------------------------------------
        function orientation = fixOrientation(this, frames)
            % What we need to do from the ellipse orientation is to figure out the head direction
            %diffi = find(abs(this.orient(frames,:) - this.direction(frames, :)) > pi/2);
            
            % The general approach is to use the presence of a tail or motion vector to disambiguate the 
            % orientation of the animal
            
            %for frames with visible tail
            tp = [this.tailVisible];
            tp = logical(tp);
            %tp(:) = 0; %let's avoid some issues for now
            %tp = (tp ==1); %make it into a logical array
            tp = tp(frames);
            if(sum(tp))
                ntf = length(frames(tp));
                ftp = frames(tp);
                tailCOM = zeros(ntf,2);
                for ii = 1:ntf
                    taili = this.tailblob(ftp(ii));
                    tailCOM(ii,:) = this.areas(ftp(ii), taili).Centroid;
                end
                tailV = tailCOM - this.bodyCOM(ftp,:); %vector towards the tail
                [theta, rho] = cart2pol(tailV(:,1), tailV(:,2));
                %we want these to be close, the opposite of the tail vector and the orientation
                od = circularDistance(this.bodyOrient(ftp), theta+pi); 
                dif = abs(od) > pi/2;  % these estimates are way off, so rotate them 180deg
                this.bodyOrient(ftp(dif)) = mod(this.bodyOrient(ftp(dif)) + pi, 2*pi);
                switchedframes = ftp(dif);
                for ii=1:length(switchedframes) %propogate that change
                    startf = switchedframes(ii);
                    frame = min(startf+1, this.nFrames);
                    odf = abs(circularDistance(this.bodyOrient(startf), this.bodyOrient(frame)));
                    while(odf > pi/2)
                        this.bodyOrient(frame) = mod(this.bodyOrient(frame) + pi, 2*pi); %rotate 180deg
                        prev = frame;
                        frame = min(frame+1, this.nFrames);
                        odf = abs(circularDistance(this.bodyOrient(prev), this.bodyOrient(frame)));
                    end
                    % now propogate change in the opposite direction
                    startf = switchedframes(ii);
                    frame = max(1, startf-1);
                    odf = abs(circularDistance(this.bodyOrient(startf), this.bodyOrient(frame)));
                    while(odf > pi/2)
                        this.bodyOrient(frame) = mod(this.bodyOrient(frame) + pi, 2*pi); %rotate 180deg
                        prev = frame;
                        frame = max(frame-1, 1);
                          odf = abs(circularDistance(this.bodyOrient(prev), this.bodyOrient(frame)));
                    end
                end

                vel_thresh = nanmax(this.vel(2:end)/2);
                vel_thresh = 200;
                tp = this.vel(2:end) > vel_thresh;  % %threshold in px/sec
                tp = logical([0; tp(:)]);
                tp = tp(frames);
                ntf = length(frames(tp));
                ftp = frames(tp);
                vel_v = this.direction(ftp);
                od = circularDistance(this.bodyOrient(ftp), vel_v); %we want these to be close
                dif = abs(od) > pi/2;  % these estimates are way off, so rotate them 180deg
                this.bodyOrient(ftp(dif)) = mod(this.bodyOrient(ftp(dif)) + pi, 2*pi);
                switchedframes = ftp(dif);
                for ii=1:length(switchedframes) %propogate that change
                    startf = switchedframes(ii);
                    frame = min(startf+1, this.nFrames);
                    odf = abs(circularDistance(this.bodyOrient(startf), this.bodyOrient(frame)));
                    while(odf > pi/2)
                        this.bodyOrient(frame) = mod(this.bodyOrient(frame) + pi, 2*pi); %rotate 180deg
                        prev = frame;
                        frame = min(frame+1, this.nFrames);
                        odf = abs(circularDistance(this.bodyOrient(prev), this.bodyOrient(frame)));
                    end
                    % now propogate change in the opposite direction
                    startf = switchedframes(ii);
                    frame = max(1, startf-1);
                    odf = abs(circularDistance(this.bodyOrient(startf), this.bodyOrient(frame)));
                    while(odf > pi/2)
                        this.bodyOrient(frame) = mod(this.bodyOrient(frame) + pi, 2*pi); %rotate 180deg
                        prev = frame;
                        frame = max(frame-1, 1);
                          odf = abs(circularDistance(this.bodyOrient(prev), this.bodyOrient(frame)));
                    end
                end
            end
            orientation = this.bodyOrient(frames);
        end
        
        % ------------------------------------------------------------------------------------------------------
        function frames = timesToFrames(this, time_range)
            % This is a utility function that returns a column vector of frame numbers for a set of time ranges,
            % specified by an n by 2 matrix of times
            if isempty(time_range) %all times
                time_range = [this.times(1) this.times(end)];
            end
            frames = [];
            for ii = 1:size(time_range,1) %the number of rows in time_range
                fi = find(this.times >= time_range(ii,1) & this.times <= time_range(ii,2));
                frames = cat(1, frames, fi(:));
            end
        end
        
        % ------------------------------------------------------------------------------------------------------
        function trimFields(this)
            % this function shortens the properties fields to nFrames
            this.COM = this.COM(1:this.nFrames,:,:);
            this.bodyCOM = this.bodyCOM(1:this.nFrames,:,:);
            this.noseblob = this.noseblob(1:this.nFrames,:);
            this.tailblob = this.tailblob(1:this.nFrames,:);
            this.tailVisible = this.tailVisible(1:this.nFrames);
            this.bodyOrient= this.bodyOrient(1:this.nFrames,:);
            this.orient= this.orient(1:this.nFrames,:,:);
            this.vel = this.vel(1:this.nFrames,:);
            this.direction= this.direction(1:this.nFrames,:);
            this.times = this.times(1:this.nFrames);
            this.areas = this.areas(1:this.nFrames,:);
            this.nblobs = this.nblobs(1:this.nFrames);
            
        end
        
        % ------------------------------------------------------------------------------------------------------
        function trimMovie(this, frames)
            % this function shortens the tracked segment and all object data to only those
            % frames specified. DO NOT USE TO SUBSAMPLE FRAMES.  WILL BREAK THINGS
            frames = intersect(1:this.nFrames, frames); %just so that some bad call doesn't f things up
            this.nFrames = length(frames);
            this.frameRange = this.frameRange(frames);
            
            this.COM = this.COM(frames,:,:);
            this.bodyCOM = this.bodyCOM(frames,:);
            this.noseblob = this.noseblob(frames,:);
            this.nosePos = this.nosePos(frames,:);
            this.tailblob = this.tailblob(frames,:);
            this.tailVisible = this.tailVisible(frames);
            this.orient= this.orient(frames,:,:);
            this.bodyOrient = this.bodyOrient(frames,:);
            this.vel = this.vel(frames,:);
            this.direction= this.direction(frames,:);
            this.times = this.times(frames);
            this.areas = this.areas(frames,:);
            this.nblobs = this.nblobs(frames);
            
        end
        
        
        % ------------------------------------------------------------------------------------------------------
        function removeStaticObjects(this)
            % function removeStaticObjects(this)
            % 
            % Goes through the areas field of detected objects and removes
            % those that are similarly present for a significant number of
            % frames.
            limit = round(this.frameRate * 3); %X sec; number of frames necessary to be eliminated as static
            dist_thresh = 1; %What do we call the same position, along each axis 
            area_thresh = 15; %the amount that an area can change frame to frame to be counted as the same
            abs_area_thresh = 30; %the upper bound for the size for deleted objects
            max_block_size = 2^11; %the maximum number of frames processed at once - need to break up for long segments
                                   %for memory reasons
            nblocks = ceil(this.nFrames/max_block_size);
            this.blobsToDelete = []; di = 0;
            for jj = 1:nblocks
                if jj == nblocks
                   sz = this.nFrames - (nblocks-1)*max_block_size;
                else
                   sz = max_block_size;
                end
                bi = (jj-1)*max_block_size+1;
                frames = (bi:(bi+sz-1))';
                framem = repmat(frames, 1, this.maxBlobs); %this is the frame number for each element
                % First of all, go through and find the common positions
                all_pos = [this.areas(frames,:).Centroid];
                all_pos = reshape(all_pos, 2, sz, this.maxBlobs); 
                %all_pos = permute(all_pos, [2 3 1]);
                %all_pos = permute(all_pos, [1 3 2]);
                all_pos = reshape(all_pos, 2, []);
                all_pos = all_pos';%makes it an blob*frame x 2 matrix, in a all frames blob1, then all frames blob2, etc order
                %all_pos = reshape(all_pos, 2, []);
                %all_pos = reshape(all_pos, [], 2); 
                nz = logical(all_pos(:,1)) | logical(all_pos(:,2)); %logical for indexing sample/blob
                nzi = find(nz);
                nz_pos = [all_pos(nz,1) all_pos(nz,2)];
                distm = ipdm(single(nz_pos));
                distm = distm + (dist_thresh+1)*eye(size(distm,1)); %makes the diagonal fall above thresh so that they aren't zero.
                dist_sel = (distm <= dist_thresh); 
                dist_rc = find(sum(dist_sel,2) > limit);
                % let's also find an area difference to make sure they are the same blob
                all_area = [this.areas(frames,:).Area];
                nz_area = all_area(nz);
                area_diff = bsxfun(@minus, nz_area, nz_area');
                area_diff = area_diff + eye(size(area_diff,1))*area_thresh; %again, offset the same area comparison
                area_sel = (abs(area_diff) < area_thresh);
                area_rc = find(sum(area_sel,2) > limit);
                %absolute area
                abs_area_rc = find(nz_area < abs_area_thresh);

                %Now, go through and find rows that are too populated
                sel = dist_sel & area_sel;
                rowcounts = sum(sel,2);

                static_rows = find(rowcounts > limit);
                static_rows2 = intersect(dist_rc, area_rc);
                static_rows2 = intersect(static_rows2, abs_area_rc);
                % now get the frames that we want to delete
                %origi = frames(nzi(static_rows2));
                for i=1:length(static_rows2)
                    origi = nzi(static_rows2(i));
                    framei = framem(origi);
                    blobi = floor((origi-1)/sz) + 1;
                    %framei = floor(origi/this.maxBlobs);
                    %blobi = origi - (framei*this.maxBlobs) + 1;
                    %blobi = floor((origi-1)/this.nFrames) + 1;
                    %framei = origi - ((blobi-1)*this.nFrames);
                    di = di+1;
                    this.blobsToDelete(di,:) = [framei, blobi];
                    %this.deleteArea(framei, blobi);
                end
            end
            for jj = 1:size(this.blobsToDelete,1)
                this.deleteArea(this.blobsToDelete(jj,1), this.blobsToDelete(jj,2));
            end
        end
        % =-------------------------------------------------------------------------------------------
        function deleteArea(this, framei, blobi)
           % function deleteArea(this, framei, blobi)
           % 
           % Just deletes the blob by saving an empty structre in its place
           props = struct('Area',0, 'Centroid',[0 0],'BoundingBox', [0 0 0 0], 'MajorAxisLength',0,...
                'MinorAxisLength',0, 'Orientation',0,'Extrema',[0 0 0 0], 'PixelIdxList',[], 'PixelList',[],'Perimeter', 0);
           this.areas(framei,blobi) = props; 
           this.nblobs(framei) = this.nblobs(framei)-1;
        end
        % ------------------------------------------------------------------------------------------------------
        function clearCalcData(this)
         % function clearCalcData(this)
         % 
         % This is a function that initializes or reintializes the computed data
            this.COM = NaN*zeros(this.nFrames, 2, this.maxBlobs);
            this.orient = NaN*zeros(this.nFrames, this.maxBlobs);
            this.vel = NaN*zeros(this.nFrames, this.maxBlobs);
            %this.noseVel = NaN*zeros(this.nFrames,1);
            %this.bodyVel = NaN*zeros(this.nFrames,1);
            this.direction = NaN*zeros(this.nFrames, this.maxBlobs);
            this.tailVisible = zeros(this.nFrames, 1);
            this.nosePos = NaN*zeros(this.nFrames, 2); %nose position
            this.nblobs = zeros(this.nFrames, 1);
            this.tailblob = NaN*zeros(this.nFrames,1);
            this.noseblob = NaN*zeros(this.nFrames,1);
            this.bodyOrient = NaN*zeros(this.nFrames, 1);
            this.bodyCOM = NaN*zeros(this.nFrames, 2);
            %initialize areas - must get these in the correct order too
            tempareas = struct('Area',0, 'Centroid',[0 0],'BoundingBox', [0 0 0 0], 'MajorAxisLength',0,...
                'MinorAxisLength',0, 'Orientation',0,'Extrema',[0 0 0 0], 'PixelIdxList',[], 'PixelList',[],'Perimeter', 0);
            this.areas = repmat(tempareas, this.nFrames, this.maxBlobs); % make a struct arraay length of nFrames
            this.fcLum = zeros(this.nFrames, 1); 
            this.blobID = NaN*zeros(this.nFrames, this.maxBlobs);
            this.blob_num = 1;
        end
        
        % ------------------------------------------------------------------------------------------------------
        function clearPathData(this)
         % function clearPathData(this)
         % Clears the path data
            this.paths = [];
        end
        % ------------------------------------------------------------------------------------------------------ 
        function detectPaths(this, time, absoluteTime, useAvgFrame)
            % Essentially just calls detectEdgesInFrame twice, once for
            % each path, and saves the results
            pb = 0;
            
            [this.paths, rew_image] = detectEdgesInFrame(this, time, absoluteTime, useAvgFrame);
            [this.paths(2), distract_image] = detectEdgesInFrame(this, time, absoluteTime, useAvgFrame);
            
            if (pb)
                order = [2 1 3];
                r = this.avgFrame; g = this.avgFrame; b = this.avgFrame;
                r(distract_image) = 255; g(distract_image) = 0; b(distract_image) = 0;
                r(rew_image) = 0; g(rew_image) = 255; b(rew_image) = 0;
                colorim = cat(3, r, g, b);
                figure;
                imshow(colorim);
                hold on;
            end
            
        end
         % ------------------------------------------------------------------------------------------------------ 
        function detectRefinePaths(this,time,absoluteTime,useAvgFrame)
            this.detectPaths(time, absoluteTime, useAvgFrame);
            this.refinePaths(1);
            this.refinePaths(2);
        end
        % ------------------------------------------------------------------------------------------------------
        function [e, eimage] = detectEdgesInFrame(this, time, absoluteTime, useAvgFrame)
            % function e = detectEdgesInFrame(this, time, absoluteTime)
            %
            % Detects the edges within an movie frame. This is the MouseTracker object, time is the time during the
            % movie and absoluteTime is a boolean flag to use movie relative time or absolute time (useful for detecting
            % edges in parts of the movie that aren't used for tracking).
            
            EDGE_LEN_THRESH = 20;
            disk_size = 20;
            pb = 0;
            % choose the image to use, then adjust to maximize contrast
            if (isempty(useAvgFrame))
                useAveFrame = 0;
            end
            if ~useAvgFrame
                if ~absoluteTime
                    f = this.timesToFrames([time time+1]);
                else % use the absolute time of the video
                    f = time*this.frameRate+1;
                    f = round(f)-this.frameRange(1);
                    f = [f f+1];
                end
                vid_struct = this.readFrames(f);
                gf = vid_struct.frames(1).cdata;
            else
                gf = this.avgFrame;
            end
            mingf = double(min(gf(:))); maxgf = double(max(gf(:)));
            gf = imadjust(gf, [mingf/255; maxgf/255], [0; 1]);
            
            % get user input about which lines they want to be marked
            fh = figure; imshow(gf);
            [cx, cy] = ginput; %get 2 points from user interaction
            close(fh);
            
            [ei, thresh] = edge(gf, 'canny'); %first detect edges
            %thresh(1) = .8*thresh(2); %inefficient, to detect twice, but I want to adjust the thresh
            %ei = edge(gf, 'canny', thresh); %first detect edges
            ei = imclose(ei, strel('square', 3)); %morphological close, fills in small (1 px) gaps
            props = {'Area', 'PixelIdxList', 'PixelList'};
            e = regionprops(ei, props); %gets the binary regions defined by the edges
            % eliminate the short edges
            ea = [e.Area];
            long = ea >= EDGE_LEN_THRESH; 
            e = e(long);
            ea = ea(long);
            [~, order] = sort(ea); e = e(order); %sort by area
            
            %find areas where we've clicked
            p = round([cx(:), cy(:)]); %clicked points
            match = zeros(length(e),1); %areas that have been matched
            for ii=1:size(p,1) 
                for jj = 1:length(e)
                    linepx = e(jj).PixelList;
                    dists = sqrt((linepx(:,1) - p(ii,1)).^2 + (linepx(:,2) - p(ii,2)).^2); %distances to every point in area
                    if min(dists) <= 10 %then we've clicked very near to some edge area
                        match(jj) =1;
                    end
                end
            end
            e = e(match == 1); % only keep areas that have been matched
            e = mergeAreas(e); % want this user interaction to result in a single area
            % This is the returned binary image
            eimage = false(size(ei));
            eimage(e.PixelIdxList) = 1;
%             eimage = imfill(eimage, 'holes'); %We want to fill in the holes in the detected path
%             e = regionprops(eimage,props); %now just redetect.
%             e = mergeAreas(e);
            
            if (pb)
                hold on;
                overlay = gf;
                neg = gf;
                overlay(e.PixelIdxList) = 255;
                neg(e.PixelIdxList) = 0;
                colorim = cat(3, overlay, neg, neg);
                imshow(colorim);
            end
        end
        % ------------------------------------------------------------------------------------------------------
        function pathIm = plotPaths(this)
            %returns a binary image (uint8 format) of the detected paths
            
            pathIm = zeros(this.height, this.width, 'uint8');
            for ii = 1:length(this.paths)
                path = this.paths(ii);
                pathIm(path.PixelIdxList) = 1;
                pathIm = imfill(pathIm,'holes');
            end
        end
        % ------------------------------------------------------------------------------------------------------
        function pathIm = plotPathsOnBg(this, pathNums)
            if (~exist('pathNums','var'))  
                pathNums = 1:length(this.paths);
            elseif(~isempty(pathNums)) %if we're given a set of paths, use that, otherwise plot all of them
                pathNums = intersect(1:length(this.paths), pathNums); %make sure we don't try to plot anything not there
            else
                pathNums = 1:length(this.paths);
            end
            pathIm = cat(3, this.avgFrame, this.avgFrame, this.avgFrame);
            %pathIm = zeros(this.height, this.width, 3);
            color_order = [2 1 3];
            for ii = 1:length(pathNums)
                %set color layer to 255
                pi = pathNums(ii);
                c = color_order(mod(pi-1, 3)+1);
                layer = pathIm(:,:,c);
                layer(this.paths(pi).PixelIdxList) = 255;
                pathIm(:,:,c) = layer;
                oc = find(1:3 ~= c);
                %set the other layers to 0
                layer = pathIm(:,:,oc(1));
                layer(this.paths(pi).PixelIdxList) = 0;
                pathIm(:,:,oc(1)) = layer;
                layer = pathIm(:,:,oc(2));
                layer(this.paths(pi).PixelIdxList) = 0;
                pathIm(:,:,oc(2)) = layer;
            end
        end
        % ------------------------------------------------------------------------------------------------------
        function refinePaths(this, pathNum)
            % function removePaths(this)
            %
            % This pops up the background with paths plotted on it (the path specified)
            % so that the user can remove portions of it.
            props = {'Area', 'PixelIdxList', 'PixelList'};
            fh = figure;
            set(fh, 'WindowKeyPressFcn', @this.exitMovieLoop);
            this.exitMovie = 0;
            while(~this.exitMovie)
                imshow(this.plotPathsOnBg(pathNum));
                del_poly = impoly(gca);
                del_roi = del_poly.createMask();
                del_i = find(del_roi);
                path_i = this.paths(pathNum).PixelIdxList;
                [keep, ki] = setdiff(path_i, del_i);
                this.paths(pathNum).PixelList = this.paths(pathNum).PixelList(ki,:);
                this.paths(pathNum).PixelIdxList = keep;
                this.paths(pathNum).Area = length(keep);   
            end 
            eimage = false(this.height, this.width);
            eimage(this.paths(pathNum).PixelIdxList) = 1;
            eimage = imfill(eimage, 'holes'); %We want to fill in the holes in the detected path
            e = regionprops(eimage,props); %now just redetect.
            this.paths(pathNum) = mergeAreas(e);
            imshow(this.plotPathsOnBg([]));
        end
        % ------------------------------------------------------------------------------------------------------
        function [noseDist, vects, ret_frames] = noseDistToTrail(this, frames, pathNum, varargin)
        % function noseDist = noseDistToTrail(this, frames, pathNum, varargin)
        % 
        % returns the distances for the nose position to the closest point of trail
        % for each of the frames specified.
        % frames - the frames to be considered
        % pathNum - path number, generally 1 for rewarded, 2 for unrewarded
        % varargin - 1) the threshold distance for which to return the frames where the distance is within.
            
            if (length(this.paths) >= pathNum)
                nosePos = this.nosePos(frames,:);
                trailPos = this.paths(pathNum).PixelList;
                distm = ipdm(single(nosePos), single(trailPos));
                [noseDist, mini] = nanmin(distm, [], 2);
                closestTrailP = trailPos(mini, :);
                vects = closestTrailP - nosePos;
            else % if there are no paths return zeros, but if there are return a mock path result
                if isempty(this.paths)
                    noseDist = NaN*zeros(length(frames), 1);
                    vects = NaN*zeros(length(frames), 2);
                    disp('There must be a detected path to calculate distances');
                else %make a fake trail
                    nosePos = this.nosePos(frames,:);
                    fakeTrail = this.makeMirrorPath(this.paths(1));
                    trailPos = fakeTrail.PixelList;
                    distm = ipdm(single(nosePos), single(trailPos));
                    [noseDist, mini] = nanmin(distm, [], 2);
                    closestTrailP = trailPos(mini, :);
                    vects = closestTrailP - nosePos;
                end
            end
            if ~isempty(varargin)
                thresh_dist = varargin{1}; %varargin 1 is the threshold distance to only return those distances under that
                seli = noseDist < thresh_dist;
                noseDist = noseDist(seli);
                ret_frames = frames(seli);
                vects = vects(seli,:);
            else
                ret_frames = frames;
            end
        end
        % ------------------------------------------------------------------------------------------------------
        function vects = noseVectorToTrail(this, frames, pathNum)
        %function vects = noseVectorToTrail(this, frames, pathNum)   
        %
        % returns the vectors for each frame. This function computes this
        % vector separately for each frame, so it may be slow if used for
        % large numbers of frames at once.  Thus it is recommended that the
        % frames be selected prior to calling this method rather than
        % performing this on a whole movie worth of frames.
            if (length(this.paths) >= pathNum)
                nosePos = this.nosePos(frames,:);
                trailPos = this.paths(pathNum).PixelList;
                for ii = 1:size(nosePos,1)
                    distm = ipdm(nosePos, trailPos);
                    nose
                end
            else
                vects = NaN*zeros(length(frames), 2);
                disp('There must be a detected path to calculate distances');
            end
            
        
        end
        % ---------------------------------------------------------------------------------------------------
        function headingVect = headingFromBodyToNose(this, frames)
        %function headingVect = headingFromBodyToNose(this, frames)   
        %
            headingVect = this.nosePos(frames,:) -  this.bodyCOM(frames,:);
            theta = cart2pol(headingVect(:,1), headingVect(:,2));
            %this.bodyOrient(frames) = theta;
            headingVect = theta;
        end
        
        % ------------------------------------------------------------------------------------------------------
        function dists = orthogonalDistFromTrail(this, frames, pathNum)
        % function dists = orthogonalDistFromTrail(this, frames)
        %
        % Let's put several pieces together
            [noseDist, vects, ~] = this.noseDistToTrail(frames, pathNum);
            orthoTheta = mod(cart2pol(vects(:,1), vects(:,2)),2*pi);
            headingTheta = mod(this.headingFromBodyToNose(frames),2*pi);
            rotations = rotationDirection(headingTheta, orthoTheta);
            %rotations = mod(headingTheta - orthoTheta, 2*pi);
            %negi = rotations > pi;
            dists = noseDist;
            dists = dists .* rotations; %gives it a sign
            %dists(negi) = -1*dists(negi); %switch the sign of some
            
        end
        
        % ------------------------------------------------------------------------------------------------------
        function meanDist = meanOrthogonalDistFromTrail(this, frames, pathNum)
        % function dists = meanOrthogonalDistFromTrail(this, frames, pathNum, threshDist)
        %
        %
            dists = this.orthogonalDistFromTrail(frames, pathNum);
            meanDist = nanmean(dists);
        
        end
        
        % ------------------------------------------------------------------------------------------------------
        function dists = orthogonalDistFromTrailPerSection(this, frames, pathNum, threshDist)
        % function dists = meanOrthogonalDistFromTrail(this, frames, pathNum, threshDist)
        %
        %
            followingFrames = this.getFollowingSegments(frames, pathNum, threshDist);
            dists = zeros(size(followingFrames,1), 1);
            for ii = 1:size(followingFrames, 1)
                dists(ii) = this.meanOrthogonalDistFromTrail(followingFrames(ii,1):followingFrames(ii,2), pathNum);
            end
        end
        
        % ------------------------------------------------------------------------------------------------------
        function plotVectsToTrail(this, frames, pathNum)
            [noseDist, vects, ~] = this.noseDistToTrail(frames, pathNum);
            np = this.nosePos(frames,:);
            
            figure;
            imshow(this.plotPathsOnBg()); %make the background
            hold on;
            quiver(np(:,1), np(:,2), vects(:,1), vects(:,2), 0);
        end
        
            
        % ------------------------------------------------------------------------------------------------------
        function pTime = propTimeOnTrail(this, frames, trailNum, threshDist)
        %function pTime = propTimeOnTrail(this, frames, trailnum)    
        %
        % Return the percentage of the time period given that the mouse 
        % was within the threshold distance from the trail.
            %traillNum = 1; % Eventually we need 2+ trails
            if nargin < 4 %default distance
                threshDist = 10;
            end
            if isempty(frames)
                frames = 1:this.nFrames;
            end
            pTime = 0; %default return value
            dists = this.noseDistToTrail(frames, trailNum);
            nn = ~isnan(dists);
            nnframes = frames(nn);
            dists = dists(nn);
            numFrames = length(nnframes);
            closeFrames = sum(dists < threshDist);
            pTime = closeFrames/numFrames;
        end
        % ---------------------------------------------------------------------------------------------------
        function followingFrames = getFollowingSegments(this, frames, trailNum, threshDist)
        % function followingFrames = getFollowingSegments(this, frames, trailNum, threshDist)
        %
        % Returns the frame segments for which the mouse is following the trail
        
            if isempty(frames)
                frames = 1:this.nFrames;
            end
            %traillNum = 1; % Eventually we need 2+ trails
            if nargin < 4 %default distance
                threshDist = 10;
            end
            distTraveled = [];
            frameNums = [];
            
            dists = this.noseDistToTrail(frames, trailNum);
            nn = ~isnan(dists);
            nnframes = frames(nn);
            dists = dists(nn);
            numFrames = length(nnframes);
            trackingInds = find(dists < threshDist); %close to trail
            trackJumps = diff(trackingInds); %identify contiguous and non are 
            skips = find(trackJumps > 1); % noncontiguous
            if ~isempty(trackingInds)
                seq_end = trackingInds(skips); %these are the indices of segments of following
                seq_start = [trackingInds(1) trackingInds(skips(1:end-1)+1)']';
                frame_start = nnframes(seq_start); %now get the actual frame numbers
                frame_end = nnframes(seq_end);
                if (isempty(frame_end)) frame_end = nnframes(trackingInds(end)); end
                frameNums = [frame_start(:) frame_end(:)];
            end
            followingFrames = frameNums;
        end
        % ---------------------------------------------------------------------------------------------------
        function [distTraveled_ret, frameNums_ret] = distanceOnTrail(this, frames, trailNum, threshDist)
        % function [dists, frames] = distanceOnTrail(this, frames, threshDist) 
        %
        % Returns the distances travelled along the trail, along with the
        % frame numbers for the following onset and offset.
        distTraveled = [];
        frameNums = getFollowingSegments(this, frames, trailNum, threshDist);
        for ii=1:size(frameNums,1)
            range = frameNums(ii,1):frameNums(ii,2);
            pos = this.nosePos(range,:);
            pos = reshape(pos(~isnan(pos)),[],2); %this is to deal with the
            % fact that there may be nan position frames in the middle.
            point_diffs = diff(pos,1,1);
            if isempty(point_diffs)
                distTraveled(ii) = 0;
            else
                distTraveled(ii) = nansum(sqrt(nansum(point_diffs.^2,2)));
            end
        end
        if isempty(distTraveled)
            distTraveled = 0;
            frameNums = [1 1];
        end
        %distTraveled_ret(ii) = distTraveled(ii);
        distTraveled_ret = distTraveled;
        frameNums_ret = frameNums;
        end
        
        % ------------------------------------------------------------------------------------------------------
        function turning_total = totalTurning(this, frames)
            measure_name = 'direction';
            measure = this.(measure_name);
            turning_vect = diff(measure);
            framei = frames > 1;
            frames = frames(framei);
            turning_vect = turning_vect(frames-1);
            turning_total = nansum(turning_vect);
        end
        % ------------------------------------------------------------------------------------------------------
        function showFrame(this, framei, movieType, dispCrop, varargin)
            % function showFrame(this, framei, useBinMovie)
            %
            % plots a frame of the movie,
            if isempty(dispCrop) dispCrop = [1 1 this.width this.height]; end
            if (nargin >= 5)
                logIm = varargin{1};
            else
                logIm = false(this.height, this.width);
            end
            dbg = 0;
            if strcmp(movieType, 'bin')
                bf = zeros(this.height, this.width, 'uint8');
                for jj=1:size(this.areas,2);
                    on = this.areas(framei,jj).PixelIdxList;
                    bf(on) = 1;
                end
                %also, need to highlight any areas to delete
                if dbg
                    fi = find(this.blobsToDelete(:,1) == framei);
                    for jj = 1:length(fi)
                        on = this.areas(this.blobsToDelete(fi(jj),1),this.blobsToDelete(fi(jj),2)).PixelIdxList;
                        bf(on) = bf(on)+jj;
                    end
                end
                pathIm = this.plotPaths()*4;
                bf = bf+pathIm + uint8(logIm)*3;
                imshow(label2rgb(bf, 'cool','k')); hold on;  
                %f = label2rgb(bf, 'cool','k');
                %pos = max(f,3); pos(logIm) = 255;
                %neg = max(f,3); f; neg(logIm) = 0;
            else
                f = this.readMovieSection(framei, movieType);
                f = this.plotPathsOnFrame(f, 1:2, 1);
                if (size(f,3) == 1)
                    f = this.gray2rgb(f);
                end
                if (nargin >=5) % A little overlay, hopefully
                    f(:,:,3) = 255 * uint8(logIm);
                    f(:,:,1) = 255 * uint8(logIm);
                end
                %neg = f; neg(logIm) = 0;
                %f = cat(3, neg, neg, pos);
                %f(:,:,3) = max(1, f(:,:,3)+cast(logIm,'uint8')*255); %not quite sure what I'm doing here
                %f(:,:,1) = 
                imshow(f);
                hold on;
            end
            %annotate the image
            if ~isempty(this.areas)
                hold on;
                for jj = 1:size(this.COM,3)
                    plot(this.COM(framei,1,jj), this.COM(framei,2,jj), 'r+', 'MarkerSize', 12, 'LineWidth',1);
                    ellipse(this.areas(framei,jj).MajorAxisLength/2, this.areas(framei,jj).MinorAxisLength/2, ...
                            this.orient(framei,jj), this.COM(framei,1,jj),this.COM(framei,2,jj),'r');
                    %text(this.COM(framei,1,jj)+10, this.COM(framei,2,jj), num2str(this.blobID(framei, jj)), 'Color','r', 'FontSize', 14);
                end
                line(this.bodyCOM(framei,1), this.bodyCOM(framei,2), 'Marker', 'o', 'Color', 'c','MarkerSize', 12, 'LineWidth',2);
                    %line(this.areas(framei).Extrema(:,1), this.areas(framei).Extrema(:,2), 'Marker', '.', 'Color', 'c');
                    %[u, v] = pol2cart(this.orient(framei), this.vel(framei)*.1);
                    %quiver(this.COM(framei,1), this.COM(framei,2), u,v, 'LineWidth', 2); %plots an orientation arrow
                if this.tailVisible(framei)
                    line(this.areas(framei,this.tailblob(framei)).Centroid(1), this.areas(framei,this.tailblob(framei)).Centroid(2), ...
                        'Marker', 'x', 'Color', 'm','MarkerSize', 12, 'LineWidth',2);
                end
                if (~isnan(this.noseblob(framei)))
                    line(this.areas(framei,this.noseblob(framei)).Centroid(1), this.areas(framei,this.noseblob(framei)).Centroid(2), ...
                        'Marker', '+', 'Color', 'g','MarkerSize', 8, 'LineWidth',2);
%                    [xv,yv] = pol2cart(this.orient(framei), this.vel(framei));
%                    quiver(this.areas(framei, this.noseblob(framei)).Centroid(1), this.areas(framei,this.noseblob(framei)).Centroid(2), xv, yv, 0);
                end
            end
           set(gca, 'xlim', [dispCrop(1) dispCrop(3)], 'ylim', [dispCrop(2) dispCrop(4)]);
        end
        
        
        % ------------------------------------------------------------------------------------------------------
        %           
        function res = readFrames(this, frames, flag)
        % function vid_struct = readFrames(this, frames, flag)
        % 
        % Flag specifies if the reading is a SINGLE frame, a CONTINUOUS range, or a DISCONTINUOUS set of frames
            adjFrames = frames+double(this.frameRange(1))-1;
            if strcmp(flag, 'single')
                res = this.readerObj.read(adjFrames(1));
                res = squeeze(res(:,:,1));
            elseif strcmp(flag,'continuous')
                res = this.readerObj.read(adjFrames);
                res = squeeze(res(:,:,1,:));
            elseif strcmp(flag, 'discontinuous') %due to limitations in videoReader class, have to loop for this
                totalFrames = sum(diff(frames')'+1);
                res = zeros(this.nativeHeight, this.nativeWidth, totalFrames); %4D matrix
                res_ind=0;
                for ii = 1:size(adjFrames,1) %loop that reads
                    ind_n = adjFrames(ii,2)-adjFrames(ii,1)+1;
                    inds = (1:ind_n) + res_ind;
                    res_ind = inds(end);
                    fi = adjFrames(ii,:);
                    if ind_n == 1
                        fi = adjFrames(ii,1);
                    end
                    temp = this.readerObj.read(fi);
                    res(:,:,inds) = squeeze(temp(:,:,1,:));
                    %res(:,:,:,inds) = this.readerObj.read(fi);
                end
            end
            res = this.applyCrop(res); %crop
        end
        
        % ------------------------------------------------------------------------------------------------------
        function movieArray = readMovieSection(this, frames, movieType, varargin)
        % function movieArray = readMovieSection(this, frames, movieType)
        % 
        % reads a section of movie specified by the frames list.
        % Inputs: frames - a list of frame ranges, formatted like [1 40; 80 100]
        %         movieType - string either 'orig' for original movie, 'diff' for the difference, or 'bin' for the
        %         binary
        %         varargin - a threshold value if you want to specify that for a binary movie 
            [frame_ints, flag] = this.listToIntervals(frames);
            rawArray = this.readFrames(frame_ints, flag);
            movieArray  = this.processFrameArray(rawArray, 1:length(frames), this.avgFrame, movieType, varargin);
        end
        
        % ------------------------------------------------------------------------------------------------------
        function movieArray = applyCrop(this, movieArray)
            % The crop is the number of px to cut in from the edge in this format [left, right, top, bottom]
            left = (this.crop(1)+1);
            horiz = (1:this.width) + left - 1;
            top = (this.crop(3)+1);
            vert = (1:this.height) + top - 1;
            movieArray = movieArray(vert, horiz, :,:);
            
        end
        
        % --------------------------------------------------
        function findMovieFile(this, vidFN)
            % function findMovieFile(this)
            %
            % Exists to find the movie file if the saved .mat file for the
            % tracking is moved since the VideoReader object saves the absolute
            % filename of the movie.
            
            disp(['Opening new video reader object for: ' vidFN]);
            this.readerObj = VideoReader(vidFN);
        end
        
        % ---------------------------------------------------
        function fakePath = makeMirrorPath(this,path)  
            fakePath.Area = path.Area;
            PixelList = [this.width - path.PixelList(:,1)+1, this.height - path.PixelList(:,2)+1];
            fakePath.PixelIdxList= this.height*(PixelList(:,1)-1) + PixelList(:,2);
            fakePath.PixelList = PixelList;
        end

        % -------------------- Functions dealing with the log file (may not exist) -------------------------------
    
        % ----------------------------------------------------------------
        function readTrackingLog(this, varargin)
        % function readTrackingLog(varargin)
        % 
        % This reads in a log file detailing the tracked segments from another program.  
        % 
            % parse input, generate default
            if nargin >= 2
                filestr = varargin{1}; % first arg is 'this'
            else
                [path, fn, ~] = fileparts(this.videoFN);
                [basename, ~] = strtok(fn, '_');
                expr = [basename, '_(.*)-0000'];
                [~,~,~,~,tokenstring,~,~]= regexp(fn,expr);
                timestr = tokenstring{1}; 
                timestr = timestr{1}; %f'ing stupid cell nesting
                filestr = [path filesep basename '_trackingLog_' timestr '.txt'];
            end
            if ~exist(filestr, 'file')
                disp('Cannot find the specified log file');
                return;
            end
            log = readCCVLog(filestr);
            this.logStruct = log;
            % Copy the path structure
            if (log.nSkelPaths > 0)
                this.paths = log.skelPaths;
                this.fullPaths = log.paths;
                w = size(log.binMovie,2); h = size(log.binMovie,1);
                for ii = 1:log.nSkelPaths
                    this.paths(ii).PixelIdxList = sub2ind([h,w], this.paths(ii).PixelList(:,2), this.paths(ii).PixelList(:,1));
                    this.fullPaths(ii).PixelIdxList = sub2ind([h,w], this.fullPaths(ii).PixelList(:,2), this.fullPaths(ii).PixelList(:,1));
                end
            elseif (log.nPaths > 0)
                this.paths = log.paths;
                w = size(log.binMovie,2); h = size(log.binMovie,1);
                for ii = 1:log.paths
                    this.paths(ii).PixelIdxList = sub2ind([h,w], this.paths(ii).PixelList(:,2), this.paths(ii).PixelList(:,1));
                end
            end
        end
        
        % ---------------------------------------------------------------
        function alignTrackingLog(this)
        % Now we have to trim it based on what limits were specified
        % when the movie was loaded.  
            startFrame = this.frameRange(1)*this.logStruct.movieSubsample;
            endFrame = ( (this.frameRange(end)-this.frameRange(1)) *this.logStruct.movieSubsample )+startFrame;
            this.logStruct.frameInfo = this.logStruct.frameInfo(startFrame:endFrame);
            nBinFrames = size(this.logStruct.binMovie,3);
            this.logStruct.binMovie = this.logStruct.binMovie(:,:, startFrame:min(nBinFrames, endFrame));
        end
        
        % ----------------------------------------------------------------
        function drawTrackedAreas(this, framei, ah)
            %function drawTrackedAreas(framei)
            if isempty(ah)
                ah = gca;
            end
            for i=1:this.logStruct.frameInfo(framei).nAreas
                area = this.logStruct.frameInfo(framei).areas(i);
                hold on;
                plot(ah, area.PixelList(:,1), area.PixelList(:,2), '.-r', 'LineWidth', 2);
            end
        end
        
        % ------------------------------------------------------------------------------------------------
        function showTrackingLogOverlayMovie(this, movieType, frames, varargin)
            % This function is for the online tracking that generates a log
            % of the position of each blob. 
            % Vargin: 1) Display crop - as in constructor
            if isempty(frames)
                frames = 1:this.nFrames;
            end
            if isfield(this.logStruct, 'binMovie');
                sub = this.logStruct.movieSubsample;
                fi = ((int32(frames)-1).*sub)+1;
                varargin{end+1} = this.logStruct.binMovie(:,:,fi);
                showMovie(this, movieType, frames, varargin{:});
            else
                disp('You must call readTrackingLog first');
            end
        end
        
        % -----------------------------------------------------------------
        function populateStructsFromLog(this) %#ok<MANU>
            % Moves the information from a separate structure to the main
            % data structures in the class.  
            
        end
            
        % -----------------------------------------------------------------
        function binMovie = makeMovieFromLog(this, frames)
            % Creates a binary movie for the requested frames from the 
            % tracking log.  
            nFrames = length(frames); %#ok<PROP>
            binMovie = false(this.height, this.width, nFrames); %#ok<PROP>
            for i = 1:nFrames %#ok<PROP>
                binMovie(:,:,i) = drawFrame(this.logStruct.frameInfo(frames(i)), this.width, this.height);
            end
        end
            
    end %Methods
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STATIC METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Static, Access = private)
        function  mov_struct = convertToGray(mov_struct)
            % function that converts the movie, frame by frame into grayscale - if it's not already
            for ii=1:length(mov_struct.frames)
                mov_struct.frames(ii).cdata = rgb2gray(mov_struct.frames(ii).cdata);
            end
        end
        
        % ---------------------------------------------------------------------------
        function avg_frame = averageFrame(mov_struct, frame_range) 
            % function to that takes the specific movie structure and computes the average frame of the specified frames
            mov = mov_struct.frames(frame_range);
            new_mov = zeros([size(mov(1).cdata,1) size(mov(1).cdata,2) length(mov)], 'uint8');
            for ii=1:length(mov)
                new_mov(:,:,ii) = squeeze(mov(ii).cdata(:,:,1));
            end
            avg_frame = uint8(round(mean(new_mov,3)));
        end
        
        % ---------------------------------------------------------------------------
        function [ret_mov, avg_frame] = processFrameArray(rawArray, frame_range, subFrame, movieType, varargin)
            % returns a  movie using the frames specified in the input. The movie can appear different and
            % is specified with MOVIETYPE: 'orig' gives the original movie, 'diff' provides a version with a 
            % frame subtracted, by default the mean frame, 'bin' is a binary thresholded image. There is
            % an optional argument 'subFrame' specifying a frame to subtract from each frame in order to
            % improve the thresholding of certain objects. 
            % RAWARRAY is in the format of dim 1,2 - height,width, 3- frame#.
            % VARARGIN{1} is the threshold level(s) for the movie, leaving
            % it out gives the default
            boostContrast = 1; %flag for boosting the contrast
            new_mov = rawArray(:,:,frame_range);
            %detection settings
            thresh(1) = .1; % the threshold level
            erode_size = 5; %the size 
            
            % make a movie from the average frame to subtract
            if ~isempty(subFrame)
                avg_frame = subFrame;
            else
                avg_frame = uint8(round(mean(new_mov,3)));
            end
            avg_mov = uint8(repmat(avg_frame, [1 1 size(new_mov,3)]));
            %diff_mov = imabsdiff(new_mov, avg_mov); %this should give a nice moving blob.
            diff_mov = new_mov - avg_mov; %this should give a nice moving blob.
            if(boostContrast)
                diff_mov = increaseMovContrast(diff_mov);
            end
            %thresh = .09; %.09; %this seems to work after image normalization, .08
            thresh(2) = 1 * graythresh(diff_mov);
            thresh(1) =  .4 * thresh(2);
            %thresh(1) = .08; .11;
            if ~isempty(varargin) && ~isempty(varargin{1}) % I'm not exactly sure MATLAB is making empty cells
                thresh(1) = varargin{1}; 
            end
            if strcmp(movieType, 'bin')
                ret_mov = false(size(diff_mov));
                nFrames = size(ret_mov,3);
                %thresh = graythresh(new_mov)*.9;
                %thresh = .5*thresh; %changing the thresh a little to include more in the blob (including more tail)
                for jj = 1
                    for ii = 1:nFrames
                        %thresh = graythresh(new_mov(:,:,ii)); % threshold for each frame separately
                        %this is to get rid of the jagged edges due to something regarding movie compression, but also
                        %removes one pixel around the edge of the contiguous sections of blob
                        bw = im2bw(diff_mov(:,:,ii),thresh(jj));
                        bw = imerode(bw, strel('square',erode_size)); 
                        ret_mov(:,:,ii) = bw;
                    end
                end
            elseif strcmp(movieType, 'diff')
                ret_mov = diff_mov;
            else % show the original movie
                if(boostContrast)
                    new_mov = increaseMovContrast(new_mov);
                end
                ret_mov = new_mov;
            end
        end
        
        % ---------------------------------------------------------------------------
        function [intervals, desc] = listToIntervals(list)
        % function intervals = listToIntervals(list)
        %
        % This function breaks up a possibly discontinuous list 
        % into a set of bounded intervals, and a flag saying if 
        % it was continuous or not.  

            % regularize list
            list = sort(list);
            list = unique(list);

            diffs = diff(list);
            breaks = find(diffs > 1);
            intervals = zeros(length(breaks)+1, 2);
            intervals(1,1) = list(1);
            for i=1:(length(breaks))
                if ~isempty(breaks) %also will be only thing that executes in only loop
                    intervals(i,2) = list(breaks(i));
                    intervals(i+1, 1) = list(breaks(i)+1);
                end
            end
            intervals(length(breaks)+1, 2) = list(end);
            if (size(intervals) > 1)
                desc = 'discontinuous';
            else
                desc = 'continuous';
            end
        end
        
        % -------------------------------------------------
        function movieType = parseMovieType(mtype_str)
            
            if strcmpi(mtype_str, 'orig')
                movieType = 0;
            elseif strcmpi(mtype_str, 'diff')
                movieType = 1;
            elseif strcmpi(mtype_str, 'bin')
                movieType = 2;
            else
                movieType = -1;
            end
        end
        
        
                
        
    end % methods - private
end %classdef