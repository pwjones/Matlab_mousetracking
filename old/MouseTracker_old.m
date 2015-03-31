classdef MouseTracker_old < handle
    properties
        MOUSE_SCALE = 10; % minimal size of the mouse used to exclude smaller objects - radius of mask
        MIN_SIZE_THRESH = 40; %minimal size of a binary area in pixels
        videoFN = []; % filename of the video
        framesPerSeg = 100 %the # of frames read at one time
        avgSubsample = 60 %sample every X frames to make the average
        width % movie dimensions
        height
        crop
        frameRate 
        totalDuration %total length of the movie
        nFrames
        trustFrameCount %mmread sometimes fails to exactly get the number of frames
        frameRange % frame numbers loaded in reference to original file
        avgFrame % the time averaged frame
        maxBlobs = 6 % number of mice tracked
        COM % center of mass for each mouse
        orient % orientation of each mouse
        vel % velocity
        direction %movement direction
        nosePos
        times % frame times
        areas % mouse regions
        tailVisible % boolean
        paths
        boostContrast = 0; %boolean switch to boost contrast of the movie
    end
    properties (GetAccess = private)
        exitMovie = 0
    end
    
    methods
        function this = MouseTracker(varargin)
            % MouseTracker(varargin) - This is an object that loads a video and then detects a mouse in that video
            % letting you get the mouse position in various parts of it.  It does this on an request basis, then 
            % saves the tracking results to minimize memory load.  
            % 
            % Args: 1 - A filename or folder of the movie to load.  If there isn't one given, then a dialog is
            % presented.
            % 2 - frame range - a vector of frame numbers to consider as the whole movie, example 1:100
            % 3 - time range - a beginning and end time to consider as the whole movie, in sec, example [60 120].
            
            % This is all to just make sure that the file is opened correctly, and that it will take a filename,
            % directory or nothing.
            movie_folder = ''; filename = '';
            if (nargin > 0)
                filename = varargin{1};
                if exist(filename, 'dir') % a folder was specified
                    movie_folder = filename;
                    filename = '';
                elseif exist(filename, 'file') % a valid filename is specified
                    this.videoFN = filename;
                end
            end
            if isempty(filename) || ~exist(filename) %a non-complete path was specified, and needs to be chosen
                [filename, movie_folder] = uigetfile([movie_folder '*.*']);
                if filename == 0 % the user has canceled the file selection
                    return;
                end
                this.videoFN = [movie_folder filename];
            else
                this.videoFN = filename;
            end
            [base_dir, base_fn, ext] =  fileparts(this.videoFN);
            mat_fn = fullfile(base_dir, [base_fn '.mat']);
            if exist(mat_fn, 'file')
                % then load the saved file
                disp('Loading the saved object');
                load(mat_fn);
            else %initialize normally 
                % Read just a couple of frames to get an idea of video speeds, etc.
                vid_struct = mmread(this.videoFN, 1); 
                vid_struct = vid_struct(1); %sometimes returning as array.  Weird.
                this.frameRate = vid_struct.rate;
                this.width = vid_struct.width;
                this.height = vid_struct.height;

                if (nargin > 3) %there is a crop vector present
                    this.crop = varargin{4};
                    this.width = vid_struct.width - sum(this.crop(1:2));
                    this.height = vid_struct.height - sum(this.crop(3:4));
                else
                    this.crop = [ 0 0 0 0];
                end

                % figure out the range of frames to be considered
                time_range = []; this.frameRange = [];
                if (nargin > 2)
                    time_range = varargin{3};
                    t_offset = vid_struct.times(1);
                    fr = round((time_range-t_offset)*this.frameRate + 1);
                    this.frameRange = fr(1):fr(2);
                elseif nargin == 2
                    this.frameRange = varargin{2};
                end
                if isempty(this.frameRange) % then specify the whole video's range
                    this.frameRange = 1:(abs(vid_struct.nrFramesTotal)-5); %undershoot to avoid potential problems
                    if vid_struct.nrFramesTotal > 0 %the frame count will be negative if there is uncertainty
                        this.trustFrameCount = 1;
                    else
                        this.trustFrameCount = 0;
                    end
                end
                % read in enough frames to create the average frame
                %avgRange = this.frameRange(1):this.avgSubsample:this.frameRange(end);
                disp('MouseTracker: reading movie to make average frame');
                %vid_struct = mmread(this.videoFN, avgRange);     
                %vid_struct = this.convertToGray(vid_struct);
                avgRange = 1:this.avgSubsample:length(this.frameRange);
                vid_struct = this.readFrames(avgRange);
                this.avgFrame = this.averageFrame(vid_struct, 1:length(vid_struct.frames));

                % Populate the fields with empty data
                this.nFrames = length(this.frameRange);
                this.times = ((1:this.nFrames)-1)/this.frameRate; % time vector starts at 0
                %this.times = ((this.frameRange(1):this.frameRange(2))-1)/this.frameRate + t_offset;
                this.clearCalcData();
            end
        end %function MouseTracker
        
        % ------------------------------------------------------------------------------------------------------
        function [wantedCOM, nose] = mousePosition(this, time_range)
        % function [wantedCOM, nose] = mousePosition(this, time_range)    
         % Returns the position of the mouse center of mass (body) and 
         % eventually the nose for the specified time range
            wantedFrames = this.timesToFrames(time_range); 
            wantedCOM = this.COM(wantedFrames,:,:);
            if isnan(wantedCOM)
                nanFrames = wantedFrames(isnan(wantedCOM(:,1,1)));
                this = this.findMouse(nanFrames);
            end 
            % get what they wanted
            wantedCOM = this.COM(wantedFrames,:,:);
            nose = [];
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
            [this, vel] = this.computeVelocity(frames);
        end
        
        % ------------------------------------------------------------------------------------------------------
        function plotVelocity(this, frames)
            % function plotVelocity(this)
            if ~exist('frames', 'var') || isempty(frames); frames = 1:this.nFrames; end
            figure;
            sorted_vel = sort(this.vel); %sorted vector of speeds for colormapping
            pathIm = this.plotPaths()*255;
            redIm = min(this.avgFrame + pathIm, 255);
            blueIm = max(this.avgFrame - pathIm, 0);
            bgIm = cat(3, blueIm, blueIm, redIm);
            imshow(bgIm); hold on;
            xlim([0 this.width]); %fit the axes to the image
            ylim([0 this.height]);
            cm = colormap(jet(this.nFrames));
            % unfortunately, I haven't been able to figure out an easier way to do colormapping of points plotted over
            % a B&W image, without the points and the image sharing the same colormap.  So, here I'm just plotting every
            % point separately, with a color indicative of the velocity of the animal.
            for ii=1:1:length(frames)
                fi = frames(ii);
                ci = find(sorted_vel == this.vel(fi), 1, 'first');
                if ~isnan(ci)
                    line('Xdata', this.nosePos(ii,1), 'Ydata', this.nosePos(ii,2), 'Marker', '.', 'MarkerSize', 6, 'Color', 'k');
                    line('Xdata', this.COM(fi,1), 'Ydata', this.COM(fi,2), 'Marker', '.', 'MarkerSize', 10, 'Color', cm(ci,:));
                end
            end 
            %Make another figure entirely to get a color scale
            fh2 = figure;
            ah2 = axes('position', [.1 .1 .7 .8], 'visible', 'off');
            cm = colormap(jet(this.nFrames));
            axes(ah2);
            pcolor([this.vel(frames), this.vel(frames)]);
            %line('parent', ah2, 'ydata', sorted_dir, 'xdata', 1:length(sorted_dir)); 
            colorbar;
        end
        
        function plotVelocityTimes(this, time_range)
            frames = this.timesToFrames(time_range);
            this.plotVelocity(frames);
        end
        
        % ------------------------------------------------------------------------------------------------------
        function plotDirection(this, frames)
            % function plotVelocity(this)
            if ~exist('frames', 'var') || isempty(frames); frames = 1:this.nFrames; end
            fh = figure;
            ah = axes('position', [.1, .1, .7 .8]);
            cm = colormap(hsv(this.nFrames));
            sorted_dir = sort(this.direction); %sorted vector for colormapping
            imshow(this.avgFrame); hold on;
            xlim([0 this.width]);
            ylim([0 this.height]);
            % unfortunately, I haven't been able to figure out an easier way to do colormapping of points plotted over
            % a B&W image, without the points and the image sharing the same colormap.  So, here I'm just plotting every
            % point separately, with a color indicative of the velocity of the animal.
            for ii=1:length(frames)
                fi = frames(ii);
                ci = find(sorted_dir == this.direction(fi), 1, 'first');
                if ~isnan(ci)
                    line('Xdata', this.COM(fi,1), 'Ydata', this.COM(fi,2), 'Marker', '.', 'MarkerSize', 10, 'Color', cm(ci,:));
                    %line('parent', ah2, 'Xdata', this.COM(ii,1), 'Ydata', this.COM(ii,2), 'Marker', '.', 'MarkerSize', 10, 'Color', cm(ci,:));
                end
            end 
            % Make another graph with the color scale
            fh2 = figure;
            ah2 = axes('position', [.1 .1 .7 .8], 'visible', 'off');
            cm = colormap(hsv(this.nFrames));
            axes(ah2);
            pcolor([sorted_dir(:), sorted_dir(:)]);
            %line('parent', ah2, 'ydata', sorted_dir, 'xdata', 1:length(sorted_dir)); 
            colorbar;
        end
        
        function plotOrientation(this, frames)
            % function plotVelocity(this)
            if ~exist('frames', 'var') || isempty(frames); frames = 1:this.nFrames; end
            figure;
            cm = colormap(hsv(this.nFrames));
            sorted_orient = sort(this.orient); %sorted vector of speeds for colormapping
            imshow(this.avgFrame); hold on;
            xlim([0 this.width]);
            ylim([0 this.height]);
            % unfortunately, I haven't been able to figure out an easier way to do colormapping of points plotted over
            % a B&W image, without the points and the image sharing the same colormap.  So, here I'm just plotting every
            % point separately, with a color indicative of the velocity of the animal.
            for ii=1:length(frames)
                fi = frames(ii);
                ci = find(sorted_orient == this.orient(fi), 1, 'first');
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
        function writeMovie(this, filename, useBinMovie)
            % function writeMovie(this, filename)
            % Writes a movie of the tracked mouse to disk at the given location
            vidWriter = VideoWriter(filename);
            open(vidWriter);
            figure;
            for ii=1:this.nFrames
                this.showFrame(ii, useBinMovie);
                currFrame = getframe; % now to the business of writing the movie
                writeVideo(vidWriter,currFrame);
                hold off;
            end
            close(vidWriter);
        end
        
        % ------------------------------------------------------------------------------------------------------
        function showMovie(this, useBinMovie, frames)
            % function showMovie(this, useBinMovie, frameRange)
            fh = figure;
            set(fh, 'WindowKeyPressFcn', @this.exitMovieLoop);
            if ~exist('frames','var') || isempty(frames)
                frames = 1:this.nFrames;
            end
            for ii=1:length(frames)
                fi = frames(ii);
                if this.exitMovie
                    break
                end
                this.showFrame(fi, useBinMovie);
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
            [basedir, fn_base, fn_ext] = fileparts(this.videoFN);
            if ~isempty(varargin)
                % save with the given filename
                fn_base = varargin{1};
            end
            full_fn = fullfile(basedir, [fn_base '.mat']);
            save(full_fn, 'this');
            disp(['Saved object in ' full_fn]);
        end
        
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRIVATE METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods %(Access = private)
        
        function this = findMouse(this, frames)
            % function this = findMouse(this, frames)
            
            % The list of properties detected in the binary image of each frame
            props = {'BoundingBox', 'Centroid', 'Perimeter', 'Area', 'Orientation', 'MajorAxisLength', ...
                     'MinorAxisLength', 'PixelList', 'PixelIdxList', 'Extrema'};
            nSegs = ceil(length(frames)/this.framesPerSeg);
            movieDone = 0;
            for jj = 1:nSegs %divide the computation up into segments so that there is never too much in memory at once
                if(movieDone); break; end
                first = (jj-1)*this.framesPerSeg + 1; %first frame of segment
                if (jj == nSegs) % the last frame
                    last = length(frames);
                else
                    last = (jj)*this.framesPerSeg;
                end
                segFrames = frames(first:last);
                
                disp(['Finding mouse in segment ' num2str(jj)]);
                frameArray = this.readMovieSection(segFrames,1);
                if size(frameArray,3) < length(segFrames) %there was a problem reading some frames, so we need to adjust
                    disp('In MouseTracker.findMouse, frameArray has too few frames, revising expectations');
                    nFramesOrig = length(segFrames);
                    segFrames = segFrames(1:size(frameArray, 3));
                    fdiff = this.nFrames - segFrames(end); %now, consider the entire length of the frames
                    %fdiff = this.nFrames - (first-1 + length(segFrames));
                    this.nFrames = this.nFrames - fdiff;
                    frames = frames(1:this.nFrames); %just trim them off the end - this error only happens for last seg
                    this.trimFields(); %trim the object property arrays
                    movieDone = 1; %don't read movie anymore
                end
                if isempty(this.areas) %initialize areas field if it hasn't been
                    areas = regionprops(frameArray(:,:,1), props); %built-in that gives stats on binary images
                    areas.tailPoints = [];
                    areas.tailPointIdx = [];
                    areas.tailCOM = [];
                    this.areas = repmat(areas, this.nFrames, this.maxBlobs); % make a struct arraay length of nFrames
                end
                for ii = 1:length(segFrames)  %have to work 1 frame at a time, unfortunately
                    lm = bwlabel(frameArray(:,:,ii)); %label matrix
                    temp_reg = regionprops(lm, props); %gives the labels for the areas
                    %E = imerode(frameArray(:,:,ii), strel('square',this.MOUSE_SCALE));
                    %e_reg = regionprops(E, props); 
                    if ~isempty(temp_reg)
                        if (1 == 0) %debug plotting
                           bin_im = zeros(this.height, this.width); 
                           %bin_im(temp_area.PixelIdxList) = 255;
                           figure; imshow(lm, [0 max(lm(:))]); colormap(jet);
                        end
                        %select the largest area blob
                        [~, maxi] = max([temp_reg.Area]); 
                        temp_reg = temp_reg(maxi);
                        temp_reg.tailPoints = []; %need to do a temp var to add the tailPoints field
                        temp_reg.tailPointIdx = [];
                        temp_reg.tailCOM=[];
                        temp_reg.Orientation = -temp_reg.Orientation; %the returned orientation is not correct, need to correct
                        this.areas(segFrames(ii)) = temp_reg;
                        if (1 == 0) %debug plotting
                           bin_im = zeros(this.height, this.width); 
                           bin_im(temp_reg.PixelIdxList) = 255;
                           figure; imshow(bin_im);
                        end
                    else %this is a kluge - unclear what the best thing to do is if you don't detect a blob
                        subI = max(1, ii-1);
                        this.areas(segFrames(ii)) = this.areas(frames(subI)); %this will error on ii=1
                    end
                end
                clear frameArray;
            end
            positions = reshape([this.areas(frames).Centroid], 2,[])';
            this.COM(frames, :) = positions;
            this.orient(frames, :) = [this.areas(frames).Orientation]./ 180 * pi; %this is the rough estimate
            this = this.computeVelocity(frames);
            this.detectTail(frames);
            this.computeOrientation(frames); % go back and get a direction/orientation
            this.nosePos(frames,:) = this.findNose(frames); %based on orientation, get the nose position
            
        end
        
        % ------------------------------------------------------------------------------------------------------
        function [this, vel] = computeVelocity(this, frames)
            % function vel = getVelocity(this, frames)
            % Computes the velocity as the frame by frame difference in position
            
            %in order to get a velocity for the first frame (0,0) position is assumed at frame 0 
            diff_com = diff(this.COM, 1, 1); %differences in COM for each x,y position
            dists = [sqrt(sum(this.COM(1,:).^2, 2)); sqrt(sum(diff_com.^2, 2))]; %euclidian distances of each frame transition
            this.vel = dists * this.frameRate; % Velocities in px/sec
            vel = this.vel(frames);
            diff_com = [this.COM(1,:); diff_com]; % add the first position to get an equal sized vector
            this.direction = cart2pol(diff_com(:,1), diff_com(:,2)); %this is the direction of motion
        end
        
        % ------------------------------------------------------------------------------------------------------
        function [this, orientation] = computeOrientation(this, frames)
            % What we need to do from the ellipse orientation is to figure out the head direction
            %diffi = find(abs(this.orient(frames,:) - this.direction(frames, :)) > pi/2);
            
            % The general approach is to use the presence of a tail or motion vector to disambiguate the 
            % orientation of the animal
            
            %for frames with visible tail
            tp = [this.tailVisible];
            tp = (tp ==1); %make it into a logical array
            tp = tp(frames);
            ntf = length(frames(tp));
            ftp = frames(tp);
            tailCOM = zeros(ntf,2);
            for ii = 1:ntf
                tailCOM(ii,:) = this.areas(ftp(ii)).tailCOM;
            end
            tailV = tailCOM - this.COM(frames(tp),:); %vector towards the tail
            [theta, rho] = cart2pol(tailV(:,1), tailV(:,2));
            %we want these to be close, the opposite of the tail vector and the orientation
            od = circularDistance(this.orient(ftp), theta+pi); 
            dif = abs(od) > pi/2;  % these estimates are way off, so rotate them 180deg
            this.orient(ftp(dif)) = mod(this.orient(ftp(dif)) + pi, 2*pi);
            switchedframes = ftp(dif);
            for ii=1:length(switchedframes) %propogate that change
                startf = switchedframes(ii);
                frame = min(startf+1, this.nFrames);
                odf = abs(circularDistance(this.orient(startf), this.orient(frame)));
                while(odf > pi/2)
                    this.orient(frame) = mod(this.orient(frame) + pi, 2*pi); %rotate 180deg
                    prev = frame;
                    frame = min(frame+1, this.nFrames);
                    odf = abs(circularDistance(this.orient(prev), this.orient(frame)));
                end
                % now propogate change in the opposite direction
                startf = switchedframes(ii);
                frame = max(1, startf-1);
                odf = abs(circularDistance(this.orient(startf), this.orient(frame)));
                while(odf > pi/2)
                    this.orient(frame) = mod(this.orient(frame) + pi, 2*pi); %rotate 180deg
                    prev = frame;
                    frame = max(frame-1, 1);
                      odf = abs(circularDistance(this.orient(prev), this.orient(frame)));
                end
            end
            
            vel_thresh = nanmax(this.vel(2:end)/2);
            vel_thresh = 300;
            tp = this.vel(2:end) > vel_thresh;  % %threshold in px/sec
            tp = logical([0; tp]);
            tp = tp(frames);
            ntf = length(frames(tp));
            ftp = frames(tp);
            vel_v = this.direction(ftp);
            od = circularDistance(this.orient(ftp), vel_v); %we want these to be close
            dif = abs(od) > pi/2;  % these estimates are way off, so rotate them 180deg
            this.orient(ftp(dif)) = mod(this.orient(ftp(dif)) + pi, 2*pi);
            switchedframes = ftp(dif);
            for ii=1:length(switchedframes) %propogate that change
                startf = switchedframes(ii);
                frame = min(startf+1, this.nFrames);
                odf = abs(circularDistance(this.orient(startf), this.orient(frame)));
                while(odf > pi/2)
                    this.orient(frame) = mod(this.orient(frame) + pi, 2*pi); %rotate 180deg
                    prev = frame;
                    frame = min(frame+1, this.nFrames);
                    odf = abs(circularDistance(this.orient(prev), this.orient(frame)));
                end
                % now propogate change in the opposite direction
                startf = switchedframes(ii);
                frame = max(1, startf-1);
                odf = abs(circularDistance(this.orient(startf), this.orient(frame)));
                while(odf > pi/2)
                    this.orient(frame) = mod(this.orient(frame) + pi, 2*pi); %rotate 180deg
                    prev = frame;
                    frame = max(frame-1, 1);
                      odf = abs(circularDistance(this.orient(prev), this.orient(frame)));
                end
            end
            
             orientation = this.orient(frames);
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
            this.COM = this.COM(1:this.nFrames,:);
            this.orient= this.orient(1:this.nFrames);
            this.vel = this.vel(1:this.nFrames);
            this.direction= this.direction(1:this.nFrames);
            this.times = this.times(1:this.nFrames);
            this.areas = this.areas(1:this.nFrames);
        end
        
        % -----------------------------------------------------------------------------------------------------
        function detectTail(this, frames)
            % Trying to detect if we can see a tail on a frame by frame basis
            
            % reset the tail detection for the frames
            scale = this.MOUSE_SCALE - 4;
            tail_thresh = 300; %the number of pixels necessary to determine that there is a tail there
            this.tailVisible(frames) = 0;
            %bin_fig = figure('position', [10 10 100 100]); ax1 = axes('parent', bin_fig);
            %tail_fig = figure('position', [110 10 100 100]); ax2 = axes('parent', tail_fig);
            % using some morphological operations to find the tail
            for ii = 1:length(frames) 
                fi = frames(ii);
                bw = zeros(this.height, this.width, 'uint8'); 
                bw(this.areas(fi).PixelIdxList) = 1; %this is the binary im
                
                % a square works pretty well for the structuring element, a cross works less well
                se = strel('square', scale); 
                % it is MUCH (>2x) faster to pack the image before doing the morphological operations
                bwp  = bwpack(bw);
                E = imerode(bwp,se, 'ispacked', this.height); %erode the image - the tail should be gone from this.
                D = imdilate(E, se, 'ispacked', this.height); % re-dilate the eroded image to reclaim the size of elements not eliminated
                D = uint8(bwunpack(D, this.height));
                T = bw-D; %these should be the tail pixels
                %T = imdilate(T, se); %get the original area (approximately) of the tail.
                tpx = find(T); %indices of tail pixels
                
                %also, need to replace the orientation of the whole mouse with body orientation, since we've already
                %done the removal process and the tail throws the orientation off
                props = {'Orientation'};
                area = regionprops(D, props);
                this.areas(fi).Orientation = -area.Orientation; %the orientation goes in the wrong direction
                this.orient(fi) = this.areas(fi).Orientation / 180 * pi; %orient is in radians
                
                if length(tpx) >= tail_thresh
                    this.tailVisible(fi) = 1;
                    tpx_y = mod(tpx-1, this.height)+1;
                    tpx_x = floor((tpx-1)/this.height) + 1;
                    tailP = [tpx_x(:) tpx_y(:)];
                    this.areas(fi).tailPointIdx = tpx;
                    this.areas(fi).tailPoints = tailP;
                    this.areas(fi).tailCOM = mean(tailP);
                end
            end
%             for ii = 1:length(frames) 
%                 fi = frames(ii);
%                 dm = squareform(pdist(this.areas(fi).Extrema)); %distance matrix
%                 mdist = median(dm,1); % the median distance
%                 outlier = mdist > 1.6*mean(mdist);
%                 if sum(outlier) % there are points that might be tail points
%                     this.tailVisible(fi) = 1;
%                     disp(['Tail visible on frame' num2str(fi)]);
%                     this.areas(fi).tailPoints = this.areas(fi).Extrema(outlier,:);
%                 end
%          end
        end
        % ------------------------------------------------------------------------------------------------------
        function nosep = findNose(this, frames)
            % The orientations for the coordinates, since the origin is in the upper-left, start with a rightward
            % vector and go clockwise.
            nosep = zeros(length(frames), 2, this.maxBlobs);
            for ii = 1:length(frames)
                fi = frames(ii);
                morient = mod(this.orient(fi), 2*pi);
                % This is the mapping between orientation and corners as detected by REGIONPROPS
                if (morient >= 0 && morient < pi/4)
                    ei = 4;
                elseif (morient > pi/4 && morient < pi/2)
                    ei = 5;
                elseif (morient > pi/2 && morient < 3*pi/4)
                    ei = 6;
                elseif (morient > 3*pi/4 && morient < pi)
                    ei = 7;
                elseif (morient >= pi && morient < 5*pi/4)
                    ei = 8;
                elseif (morient >= 5*pi/4 && morient < 3*pi/2)
                    ei = 1;
                elseif (morient >= 3*pi/2 && morient < 7*pi/4)
                    ei = 2;
                elseif (morient >= 7*pi/4 && morient < 2*pi)
                    ei = 3;
                end
                corner = this.areas(fi).Extrema(ei,:);
                corner_vect = corner - this.COM(fi,:);
                % So, 
                nosep(ii,:) = this.areas(fi).Extrema(ei,:);
            end
        end
        % ------------------------------------------------------------------------------------------------------
        function clearCalcData(this)
         % function clearCalcData(this)
         % 
         % This is a function that initializes or reintializes the computed data
            this.COM = NaN*zeros(this.nFrames, 2, this.maxBlobs);
            this.orient = NaN*zeros(this.nFrames, this.maxBlobs);
            this.vel = NaN*zeros(this.nFrames, this.maxBlobs);
            this.direction = NaN*zeros(this.nFrames, this.maxBlobs);
            this.tailVisible = zeros(this.nFrames, this.maxBlobs);
            this.nosePos = NaN*zeros(this.nFrames, 2, this.maxBlobs); %nose position
            this.areas = [];
        end
        % ------------------------------------------------------------------------------------------------------
        function [e, eimage] = detectEdgesInFrame(this, time, absoluteTime)
            % function e = detectEdgesInFrame(this, time, absoluteTime)
            %
            % Detects the edges within an movie frame. This is the MouseTracker object, time is the time during the
            % movie and absoluteTime is a boolean flag to use movie relative time or absolute time (useful for detecting
            % edges in parts of the movie that aren't used for tracking).
            
            EDGE_LEN_THRESH = 20;
            disk_size = 20;
            
            if ~absoluteTime
                f = this.timesToFrames([time time+1]);
            else % use the absolute time of the video
                f = time*this.frameRate+1;
                f = round(f)-this.frameRange(1);
                f = [f f+1];
            end
            vid_struct = this.readFrames(f);
            gf = vid_struct.frames(1).cdata;
            %gf = 255 - gf; %invert the image
            
            % Filter out the background in order to detect sharp changes
%             sigma = 25;
%             h = fspecial('gaussian', 6*sigma, sigma);
%             background = imfilter(gf, h);
%             gf = gf - background;
            % adjust the image to take the whole range
            mingf = double(min(gf(:))); maxgf = double(max(gf(:)));
            gf = imadjust(gf, [mingf/255; maxgf/255], [0; 1]);
            %background = imopen(gf,strel('disk', disk_size)); %finds a background via opening of the grayscale image
            %gf2 = gf-background;
            
            % this step is to get rid of the bright highlights in the image.
            %medgf = median(double(gf(:))) + 30;
            %i = gf > medgf;
            %gf(i) = medgf;
            
            %ei = edge(gf2, 'canny');
            %ei = imfill(ei, 'holes');
            figure; imshow(gf);
            [cx, cy] = ginput; %get 2 points from user interaction
            
            ei = edge(gf, 'canny');
            ei = imclose(ei, strel('square', 3));
            %ei = imfill(ei);
            props = {'Area', 'PixelIdxList', 'PixelList'};
            e = regionprops(ei, props);
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
            e = e(match == 1); %only keep areas that have been matched
            
            % Make a new image from only those
            ei = zeros(size(ei));
            for ii=1:length(e)
                ei(e(ii).PixelIdxList) = 1;
            end
            eimage = (ei==1);
            eimage = imclose(eimage, strel('square', 5));
            
            this.paths = e;
            %traces = bwtraceboundary(eimage, p, 'e');
        end
        
        function pathIm = plotPaths(this)
            %returns a binary image (uint8 format) of the detected paths
            
            pathIm = zeros(this.height, this.width, 'uint8');
            for ii = 1:length(this.paths)
                path = this.paths(ii);
                pathIm(path.PixelIdxList) = 1;
            end
        end
        % ------------------------------------------------------------------------------------------------------
        function showFrame(this, ii, useBinMovie)
            if useBinMovie
                on = this.areas(ii).PixelIdxList;
                bf = zeros(this.height, this.width, 'uint8');
                bf(on) = 1;
                if this.tailVisible(ii)
                    bf(this.areas(ii).tailPointIdx) = 2;
                end
                pathIm = this.plotPaths()*4;
                bf = bf+pathIm;
                imshow(bf, [0 0 0; 1 1 1; 1 0 0; 0 0 1; 1 1 1; 1 0 0]); hold on;  
            else
                %f = mmread(this.videoFN, this.frameRange(ii));
                %f = this.readFrames(ii);
                %imshow(f.frames.cdata); 
                f = this.readMovieSection(ii,0);
                imshow(f);
                %imshow(imabsdiff(f.frames.cdata, this.avgFrame));
                hold on;
            end
            %annotate the image
            if ~isempty(this.areas) 
                plot(this.COM(ii,1), this.COM(ii,2), 'r+', 'MarkerSize', 12, 'LineWidth',1);
                ellipse(this.areas(ii).MajorAxisLength, this.areas(ii).MinorAxisLength, this.orient(ii), this.COM(ii,1),this.COM(ii,2),'r');
                line(this.areas(ii).Extrema(:,1), this.areas(ii).Extrema(:,2), 'Marker', '.', 'Color', 'c');
                [u, v] = pol2cart(this.orient(ii), this.vel(ii)*.1);
                quiver(this.COM(ii,1), this.COM(ii,2), u,v, 'LineWidth', 2); %plots an orientation arrow
                if this.tailVisible(ii)
                    line(this.areas(ii).tailCOM(1), this.areas(ii).tailCOM(2), 'Marker', 'x', 'Color', 'c','MarkerSize', 12, 'LineWidth',2);
                end
                plot(this.nosePos(ii,1), this.nosePos(ii,2), '.y', 'MarkerSize', 12);
            end
        end
        
        % ------------------------------------------------------------------------------------------------------
        function vid_struct = readFrames(this, frames)
            vid_struct = mmread(this.videoFN, frames+this.frameRange(1)-1); %this final call is where we compensate for an offset 
            vid_struct = this.convertToGray(vid_struct);
            vid_struct = this.applyCrop(vid_struct);
        end
        
        % ------------------------------------------------------------------------------------------------------
        function movieArray = readMovieSection(this, frames, binary)
            vid_struct = this.readFrames(frames);
            if binary
                movieArray = this.convertMovieToFrameArray(vid_struct, 1:length(frames), this.avgFrame, 1);
            else
                movieArray = this.convertMovieToFrameArray(vid_struct, 1:length(frames), this.avgFrame, 0);
            end
        end
        
        % ------------------------------------------------------------------------------------------------------
        function movieStruct = applyCrop(this, movieStruct)
            % The crop is the number of px to cut in from the edge in this format [left, right, top, bottom]
            left = (this.crop(1)+1);
            horiz = (1:this.width) + left - 1;
            top = (this.crop(3)+1);
            vert = (1:this.height) + top - 1;
            for ii = 1:length(movieStruct.frames)
               movieStruct.frames(ii).cdata = movieStruct.frames(ii).cdata(vert, horiz); 
            end
        end 
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STATIC METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Static, Access = private)
        function  mov_struct = convertToGray(mov_struct)
            % function that converts the movie, frame by frame into grayscale - if it's not already
            for ii=1:length(mov_struct.frames)
                mov_struct.frames(ii).cdata = rgb2gray(mov_struct.frames(ii).cdata);
            end
        end
        
        function avg_frame = averageFrame(mov_struct, frame_range) 
            % function to that takes the specific movie structure and computes the average frame of the specified frames
            mov = mov_struct.frames(frame_range);
            new_mov = zeros([size(mov(1).cdata,1) size(mov(1).cdata,2) length(mov)], 'uint8');
            for ii=1:length(mov)
                new_mov(:,:,ii) = squeeze(mov(ii).cdata(:,:,1));
            end
            avg_frame = uint8(round(mean(new_mov,3)));
        end
        
        function [ret_mov, avg_frame] = convertMovieToFrameArray(mov_struct, frame_range, subFrame, binary)
            % returns a binary thresholded movie using the frames specified in the input. Also there is
            % an optional argument 'subFrame' specifying a frame to subtract from each frame in order to
            % improve the thresholding of certain objects.
            boostContrast = 1; %flag for boosting the contrast
            if (length(mov_struct.frames) < length(frame_range))
                disp('In MouseTracker.convertMovieToBinaryArray: frame_range specified is too large for movie.  Trimming');
                frame_range = frame_range(1:length(mov_struct.frames));
            end
            mov = mov_struct.frames(frame_range);
            new_mov = zeros([size(mov(1).cdata,1) size(mov(1).cdata,2) length(mov)], 'uint8');
            for ii=1:length(mov)
                new_mov(:,:,ii) = squeeze(mov(ii).cdata(:,:,1));
            end
            % make a movie from the average frame to subtract
            if ~isempty(subFrame)
                avg_frame = subFrame;
            else
                avg_frame = uint8(round(mean(new_mov,3)));
            end
            avg_mov = repmat(avg_frame, [1 1 length(mov)]);
            new_mov = imabsdiff(new_mov, avg_mov); %this should give a nice moving blob.
            if(boostContrast)
                new_mov = increaseMovContrast(new_mov);
            end
            if binary
                ret_mov = zeros(size(new_mov), 'uint8');
                nFrames = size(ret_mov,3);
                thresh = graythresh(new_mov)*.9;
                %thresh = .5*thresh; %changing the thresh a little to include more in the blob (including more tail)
                for ii = 1:nFrames
                    %thresh = graythresh(new_mov(:,:,ii)); % threshold for each frame separately
                    %this is to get rid of the jagged edges due to something regarding movie compression, but also
                    %removes one pixel around the edge of the contiguous sections of blob
                    bw = im2bw(new_mov(:,:,ii),thresh);
                    %bwf = imerode(bw, strel('square',3)); 
                    ret_mov(:,:,ii) = bw;
                end
            else
                ret_mov = new_mov;
            end
        end
        
        
    end % methods - private
end %classdef