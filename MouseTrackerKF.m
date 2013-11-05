classdef MouseTrackerKF < MouseTracker
    properties
        default_thresh = 0;
        used_thresh = 0;
        p_mouse = .0008; %the proportion of pixels that are 
        kf = struct('nstates', {}, 'nob', {}, 's', {});
    end
    
    methods
        function this = MouseTrackerKF(varargin)
            % MouseTrackerKF(varargin) - This is an object that loads a video and then detects a mouse in that video
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
            
            % This is all to just make sure that the file is opened correctly, and that it will take a filename,
            % directory or nothing.
            
            % All of the proper initialization is done in the MouseTracker object this is inherited from
            this = this@MouseTracker(varargin{:});
%            
        end %function MouseTrackerKF
        % ------------------------------------------------------------------------------------------------------
        function initKalmanFilter(this)
        % function initKalmanFilter(this)
        %
        % Initialization of Kalman filter parameters
            this.kf = struct('nstates', {}, 'nob', {}, 's', {});
            zs = 4;
            this.kf(1).nstates = zs;
            this.kf(1).nob = 2;
            this.kf(1).s.A = eye(zs);
            this.kf(1).s.A = [1 0 1 0; 0 1 0 1; 0 0 1 0; 0 0 0 1]; 
            this.kf(1).s.z = nan*zeros(zs,1); 
            this.kf(1).s.x = nan*zeros(zs,1); 
            this.kf(1).s.H = eye(zs);
            this.kf(1).s.R = eye(zs); %measurement error covariance
            %this.kf.s.P = cov(np');
            this.kf(1).s.P = eye(zs);
            this.kf(1).s.u = zeros(zs,1);
            this.kf(1).s.B = eye(zs);
            this.kf(1).s.Q = eye(zs);

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
            imshow(this.plotPathsOnBg()); hold on;
            xlim([0 this.width]); %fit the axes to the image
            ylim([0 this.height]);
            %cm = colormap(spring(this.nFrames));
            vel_vect = this.vel(frames,1);
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
                    % So, this way we can't plot a solid line.
                    line('Parent', ah, 'Xdata', this.nosePos(fi,1), 'Ydata', this.nosePos(fi,2), ...
                        'Marker', '.','MarkerSize', 8, 'Color', c);
                end
            end 
        end
        % ------------------------------------------------------------------------------------------------------
        function plotFilterPosition(this,frames)
        % function plotFilterPosition(this,frames)
            plotblue = [.3, .6, 1]; % a nice blue color
            if ~exist('frames', 'var') || isempty(frames); frames = 1:this.nFrames; end
            fh = figure;
            ah = axes('position', [.1, .1, .7 .8]);
            bgIm = this.plotPathsOnBg();
            imshow(bgIm); hold on;
            xlim([0 this.width]);
            ylim([0 this.height]);
            filt_x = [this.kf.s.x]; %the kalman filter positions
            filt_x = filt_x(1:2,:)';
            for ii=1:length(frames)
                fi = frames(ii);
                if ~isnan(this.noseblob(fi))
                    line('Parent', ah, 'Xdata', this.areas(fi,this.noseblob(fi)).Centroid(1), 'Ydata', ...
                        this.areas(fi,this.noseblob(fi)).Centroid(2), 'Marker', '.', 'MarkerSize', 8, 'Color', plotblue);
                    line('Parent', ah, 'Xdata', filt_x(fi,1), 'Ydata', filt_x(fi,2), 'Marker', '.', 'MarkerSize', 8, 'Color', [1 0 0]);
                end
            end 
            title(this.videoFN);
        end
        % ------------------------------------------------------------------------------------------------------
        function plotNosePosition(this, frames)
            % function plotDirection(this, frames)
            plotblue = [.3, .6, 1]; % a nice blue color
            if ~exist('frames', 'var') || isempty(frames); frames = 1:this.nFrames; end
            fh = figure;
            ah = axes('position', [.1, .1, .7 .8]);
            bgIm = this.plotPathsOnBg();
            %bgIm = this.plotPaths();
            imshow(bgIm); hold on;
            xlim([0 this.width]);
            ylim([0 this.height]);
            for ii=1:length(frames)
                fi = frames(ii);
                if ~isnan(this.noseblob(fi))
                    %line('Parent', ah, 'Xdata', this.areas(fi,this.noseblob(fi)).Centroid(1), 'Ydata', ...
                    %    this.areas(fi,this.noseblob(fi)).Centroid(2), 'Marker', '.', 'MarkerSize', 8, 'Color', plotblue);
                    line('Parent', ah, 'Xdata', this.nosePos(fi,1), 'Ydata', ...
                        this.nosePos(fi,2), 'Marker', '.', 'MarkerSize', 10, 'Color', plotblue);
                    %line('Parent', ah, 'Xdata', this.nosePos(fi,1), 'Ydata', ...
                    %    this.nosePos(fi,2), 'Marker', '.', 'MarkerSize', 10, 'Color', 'w');
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
            [cm cinds] = getIndexedColors('jet', this.direction(frames), 1);
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
                    line('Xdata', this.nosePos(fi,1), 'Ydata', ...
                        this.nosePos(fi,2), 'Marker', '.', 'MarkerSize', 8, 'Color', cm(cinds(ii),:));
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
            colormap(cm);
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
        
        % ------------------------------------------------------------------------------------------------------
        function showMovie(this, movieType, frames, varargin)
            % function showMovie(this, useBinMovie, frameRange)
            if length(varargin) >=1
                dispCrop = varargin{1};
            else
                dispCrop = [];
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
                this.showFrame(fi, movieType, dispCrop);
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
        
        function mov = returnFrames(this, frames, movieType)
            % mov = returnFrames(this, frames, binary)
            % returns the frames of the movie specified, binary or greyscale
            mov = this.readMovieSection(frames, movieType);
            
        end
        
        % Setting functions for manually altering the tracking result
        function setTailPosition(this, frame, pos)
            frame = frame(1);
            centers = [this.areas(frame, :).Centroid];
            centers = reshape(centers, 2, [])';
            dist = ipdm(centers, pos);
            if nanmin(dist) < 10
                [min_dist, disti] = nanmin(dist);
                this.tailblob(frame) = disti;
            end
        end
        
        function setNosePosition(this, frame, pos)
            frame = frame(1);
            centers = [this.areas(frame, :).Centroid];
            centers = reshape(centers, 2, [])';
            dist = ipdm(centers, pos);
            if nanmin(dist) < 10
                [min_dist, disti] = nanmin(dist);
                this.noseblob(frame) = disti;
                this.nosePos(frame,:) = this.areas(frame,disti).Centroid;
            end
        end
        
    %end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRIVATE METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %methods %(Access = private)
        
        function this = findMouse(this, frames)
            % function this = findMouse(this, frames)
            
            nSegs = ceil(length(frames)/this.framesPerSeg);
            movieDone = 0; 
            
            
            for jj = 1:nSegs %divide the computation up into segments so that there is never too much in memory at once
                first = (jj-1)*this.framesPerSeg + 1; %first frame of segment
                if (jj == nSegs) % the last frame
                    last = length(frames);
                else
                    last = (jj)*this.framesPerSeg;
                end
                segFrames = frames(first:last);
                
                disp(['Finding mouse in segment ' num2str(jj)]);
                frameArray = this.readMovieSection(segFrames,'bin', this.default_thresh);
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
                detectMouseInSection(this, segFrames, frameArray);
                clear frameArray;
            end
            if trimmed %means the last set of frames came out shorter than expected from movie info, trimming things
                this.nFrames = min(newFrameCount);
                frames = frames(1:this.nFrames); %just trim them off the end - this error only happens for last seg
                this.trimFields(); %trim the object property arrays
            end
            this.computeVelocity(frames);
            %this.refineTracking(frames);    
        end
        
        % -----------------------------------------------------------------------------------------------
        function detectMouseInSection(this, segFrames, frameArray)
            
            err_thresh = 6;
            dbg = 0;
            cb = 0;
            
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
                    this.nblobs(fi) = min(this.maxBlobs, length(temp_reg)); %take the biggest N blobs
                    this.areas(segFrames(ii),1:this.nblobs(fi)) = temp_reg(1:this.nblobs(fi)); %only keep the top sized blobs
                    
                    this.detectTail(fi);
                    this.bodyCOM(fi,:) = this.computeBodyPos(fi,1); % 1) get the center of mass of animal
                    [this.nosePos(fi,:), this.noseblob(fi)] = this.findNose(fi); % 2) get the nose blob
                    s = this.kf.s(fi); % 3) compare to prediction
                    pred_x = s.x(1:2)';
                    if( fi == 1) %if first frame, make sure that things are good
                        this.kf.s(fi).z = [this.nosePos(fi,:) 0 0]';
                        pred_x = this.kf.s(fi).z(1:2)';
                    end
                    if (pdist2(pred_x, this.nosePos(fi,:)) > err_thresh) || isnan(this.nosePos(fi,1))
                        % The difference between prediction and detected position is too big, so try
                        % to correct (or the nose wasn't found)
                        if cb
                            this.correctDetection(fi, this.nosePos(fi,:), pred_x);
                        end
                    end
                    % just update the filter with the proper z in order to predict next frame
                    if fi>1 %compute frame-by-frame velocities
                        this.bodyVel(fi, :) = this.bodyCOM(fi,:) - this.bodyCOM(fi,:);
                        this.noseVel(fi, :) = this.nosePos(fi,:)-this.nosePos(fi-1,:);
                        if isnan(this.nosePos(fi-1,1))
                            vel = [0 0];
                        else
                            vel = this.noseVel(fi, :);
                        end
                    else
                        vel = [0 0];
                    end
                    this.kf.s(fi).z = [this.nosePos(fi,:) vel]';
                    % 5) make prediction for next frame
                    if fi < this.nFrames
                        this.kf.s(fi+1) = kalmanf(this.kf.s(fi));
                    end
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
                    this.areas(segFrames(ii)) = this.areas(segFrames(subI)); %this will error on ii=1
                    if fi < this.nFrames
                        this.kf.s(fi+1) = kalmanf(this.kf.s(fi)); %let's update the Kalman Filter anyway
                    end
                end
            end
        end
        % ------------------------------------------------------------------------------------------------------
        function newPos = correctDetection(this, frame, detectedPos, predictedPos)
            err_thresh = 7;
            areas = this.areas(frame,:);
            centers = [areas.Centroid];
            centers = reshape(centers, 2,[])';
            D = pdist2(predictedPos, centers); %distances from prediction to all centers
            [mind, mini] = nanmin(D);
            prevPos = [0 0];
            if frame > 1
                prevPos = this.nosePos(frame-1,:);
            end
            if isnan(this.noseblob(frame)) && mind < err_thresh %if we missed the nose blob for some reason
                newPos = centers(mini, :);
                this.nosePos(frame, :) = newPos;
                this.noseblob(frame) = mini;
                this.kf.s(frame).z = [newPos newPos-prevPos]';
                disp(sprintf('Corrected frame %i from %d, %d   to %d, %d', frame, detectedPos(1), detectedPos(2), newPos(1), newPos(2)));
            elseif mini ~= this.noseblob(frame) && mind < err_thresh %if a different blob is closer to the prediction
                % then we will call that the nose instead
                newPos = centers(mini, :);
                this.nosePos(frame, :) = newPos;
                this.noseblob(frame) = mini;
                this.kf.s(frame).z = [newPos newPos-prevPos]';
                disp(sprintf('Corrected frame %i from %d, %d   to %d, %d', frame, detectedPos(1), detectedPos(2), newPos(1), newPos(2)));
            elseif mini == this.noseblob(frame) && mind > err_thresh
                % do nothing for now
                % eventually want to rethreshold the frame
                
            end
            
        end

        % ------------------------------------------------------------------------------------------------------
        function propogateNosePosition
            
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
                fi = frames(ii);
                regions = this.areas(fi,:);
                maxo = 0;
                if (fi>1) && this.tailVisible(fi-1) %finding similarity with previously id'd tail 
                    prev_tail = this.areas(fi-1,this.tailblob(fi-1));
                    overlap = regionOverlap(prev_tail, regions, this.regprops);
                    [maxo, maxoi] = max(overlap);
                end
                if maxo > .1 %if similar, we've found another tail - this is a conservative threshold. Comparing the 
                    % overlap of id'd tail blobs to non-tail blobs, they are all down below
                    % .05...essentially zero.  This is really conservative - assuming that you get the
                    % tail in the first place.
                    this.tailVisible(fi) = 1;
                    this.tailblob(fi) = maxoi;
                else
                    % going to find the blob with the highest eccentricity to be the tail
                    ecc = [regions.MajorAxisLength]./[regions.MinorAxisLength];
                    [max_ecc, maxi] = nanmax(ecc);
                    if max_ecc>4.5 && regions(maxi).Area > 100%this is just an empirically defined threshold value
                        this.tailVisible(fi) = 1;
                        this.tailblob(fi) = maxi;
                    else
                        this.tailVisible(fi) = 0;
                        this.tailblob(fi) = NaN;
                    end
                end
            end
        end
        % ------------------------------------------------------------------------------------------------------
        function [nosep, noseblob] = findNose(this, frames)
        % function nosep = findNose(this, frames, blobi)
        %
        % So, the method for finding which blob is the nose is to take the tail and body center and take the 
        % distance along that vector for every blob centroid. The one with the largest is the nose.
            %nosep = NaN* zeros(length(frames), 2);
            nosep = this.nosePos(frames,:);
            noseblob = this.noseblob(frames);
            for ii = 1:length(frames)
                fi = frames(ii);
                regions = this.areas(fi,:);
                maxo = 0;
                if fi >1 && ~isnan(this.noseblob(fi-1)) %finding similarity with previously id'd tail 
                    if (ii == 1) prev_nb = this.noseblob(fi-1);
                    else prev_nb = noseblob(ii-1);
                    end
                    prev_nose = this.areas(fi-1,prev_nb);
                    overlap = regionOverlap(prev_nose, regions, this.regprops);
                    [maxo, maxoi] = max(overlap);
                end
                if maxo > .7 %nearly all non-nose blobs have zero overlap with a well defined nose
                    noseblob(ii) = maxoi;
                elseif (this.nblobs(fi) > 3) %only works if there are enough things to track
                    if this.tailVisible(fi) && (this.nblobs(fi) > 2)%we can only do this first method if the tail is present
                        % compute the vector from tail to body center
                        bodyVect = this.bodyCOM(fi,:) - this.areas(fi,this.tailblob(fi)).Centroid;
                        dist = [];
                        for jj=1:this.nblobs(fi)
                            areaVect = this.areas(fi,jj).Centroid - this.areas(fi,this.tailblob(fi)).Centroid;
                            dist(jj) = dot(areaVect, bodyVect);
                        end
                        [maxd, maxi] = max(dist); %the maximum distance and ind along tail-body vector of each blob  
                        noseblob(ii) = maxi; 
                        nosep(ii,:) = this.areas(fi,noseblob(ii)).Centroid; % The nose is the farthest along vector
                    end
                end
            end
            % putting in a maximum distance criterion of 50 px away from last frame 
            
            %let's trim the vectors
            %nosep = nosep(frames,:);
            %noseblob = noseblob(frames);
            
            dists = diff(nosep,1,1);
            dists = sqrt(sum(dists.^2,2));
            toofar = dists > 50;
            nosep(toofar,:) = NaN*zeros(sum(toofar), 2);
        end
        
        % ------------------------------------------------------------------------------------------------------
        function bodyCOM = computeBodyPos(this, frames, includeTail)
            %
            if isempty(includeTail) includeTail=0; end
            bodyCOM = NaN*zeros(length(frames), 2);
            for i = 1:length(frames)
                fi = frames(i);
                taili = this.tailblob(fi);
                if (isnan(taili)) taili = this.maxBlobs + 1; end %just assign it out of range for check
                temp_pos = [];
                temp_areas = this.areas(fi,1:this.nblobs(fi));
                for j=1:this.nblobs(fi)
                    this.orient(fi, j) = [this.areas(fi,j).Orientation]./ 180 * pi; %this is the rough estimate
                    positions = permute(reshape([temp_areas(j).Centroid], 2,[]), [3 1 2]);
                    this.COM(fi, :, j) = positions;
                    if (includeTail) || (j ~= taili) %exclude the tail in the body position calculation
                        %temp_pos = cat(1, temp_pos, temp_areas(j).PixelList);
                        temp_pos = cat(1, temp_pos, positions);
                        
                        %this.orient(fi, j) = [this.areas(fi,j).Orientation]./ 180 * pi; %this is the rough estimate
                    end
                end
                %bodyCOM(i, :) = mean(this.COM(fi,:,:),3);
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
            %this.vel = NaN*zeros(size(this.bodyCOM,1),size(this.bodyCOM,3));
            if frames(1) > 1 %compute velocity from the frame before the range to get all frames requested
                fr = [frames(1)-1; frames(:)];
            else
                fr = [1; frames(:)];
            end
            %pos = this.bodyCOM(fr,:);
            pos = this.nosePos(fr,:);
            diff_pos = diff(pos, 1,1);%differences in each x,y position
            vel = sqrt(sum(diff_pos.^2, 2));
            this.vel(fr(2:end),1) = vel; % vel is 1 shorter than the fr
            this.direction = cart2pol(diff_pos(:,1), diff_pos(:,2)); %this is the direction of motion
            
            %diff_com = [COM(1,:); diff_com]; % add the first position to get an equal sized vector
            %this.direction(frames) = cart2pol(diff_com(:,1), diff_com(:,2)); 
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
            this.vel = NaN*zeros(this.nFrames, 1);
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
            this.initKalmanFilter();

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
            %pathIm = cat(3, this.avgFrame, this.avgFrame, this.avgFrame);
            pathIm = zeros(this.height, this.width, 3);
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
                seq_end = [trackingInds(skips)' trackingInds(end)]'; %these are the indices of segments of following
                seq_start = [trackingInds(1) trackingInds(skips+1)']';
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
        function prop = propTrailFollowed(this, frames, trailNum, threshDist)
        % function propTrailFollowed(this, frames, trailNum, threshDist)
        %
        % This function returns the proporation of the trail pixels that the mouse came within a given
        % distance of.  Seems like a clean way of measuring the completeness of his trail exploration.
        % Algorthmically, this is a similar problem to finding the following segments except reversed.
        if isempty(frames)
            frames = 1:this.nFrames;
        end
        np = this.nosePos(frames,:);
        nn = ~isnan(np(:,1));
        np = np(nn,:);
        trailPos = this.paths(trailNum).PixelList; 
        distm = ipdm(single(np), single(trailPos));
        [trailDist, mini] = nanmin(distm, [], 1);
        explored = find(trailDist <= threshDist);
        npx = length(this.paths(trailNum).PixelIdxList);
        prop = length(explored)/npx;
        
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
        function nose_jumps = findNoseJumps(this, dist_thresh, frames)
            % reports the first frame after there is a large jump in the nose position
            if isempty(frames)
                frames = 1:this.nFrames;
            end
            dists = sqrt(sum(diff(this.nosePos(frames,:)).^2, 2));
            nose_jumps = find(dists >= dist_thresh);
        end
            
        % ------------------------------------------------------------------------------------------------------
        function showFrame(this, framei, movieType, dispCrop)
            % function showFrame(this, framei, useBinMovie)
            %
            % plots a frame of the movie,
            if isempty(dispCrop) dispCrop = [1 1 this.width this.height]; end
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
                bf = bf+pathIm;
                imshow(label2rgb(bf, 'cool','k')); hold on;
                %imshow(bf, [0 0 0; 1 1 1; 1 0 0; 0 0 1; 1 1 1; 1 0 0]); hold on;  
            else
                %f = mmread(this.videoFN, this.frameRange(framei));
                %f = this.readFrames(framei);
                %imshow(f.frames.cdata); 
                f = this.readMovieSection(framei, movieType);
                imshow(f);
                %imshow(imabsdiff(f.frames.cdata, this.avgFrame));
                hold on;
            end
            %annotate the image
            if ~isempty(this.areas)
                hold on;
                for jj = 1:size(this.COM,3)
                    plot(this.COM(framei,1,jj), this.COM(framei,2,jj), 'r+', 'MarkerSize', 12, 'LineWidth',1);
                    ellipse(this.areas(framei,jj).MajorAxisLength/2, this.areas(framei,jj).MinorAxisLength/2, ...
                            this.orient(framei,jj), this.COM(framei,1,jj),this.COM(framei,2,jj),'r');
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
                        'Marker', '+', 'Color', 'g','MarkerSize', 12, 'LineWidth',2);
                    [xv,yv] = pol2cart(this.orient(framei), this.vel(framei));
                    quiver(this.areas(framei, this.noseblob(framei)).Centroid(1), this.areas(framei,this.noseblob(framei)).Centroid(2), xv, yv, 0);
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
            movieArray = this.processFrameArray(rawArray, 1:length(frames), this.avgFrame, movieType, varargin{:});
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
        
        % ---------------------------------------------------
        function rethreshold(this, frameRange, threshold)
        % Setting the threshold differently in order to improve blob identification
            this.default_thresh = threshold;
            mArray = this.readMovieSection([frameRange(1) frameRange(end)], 'bin');
            
        end
        % ---------------------------------------------------------------------------
        function [ret_mov, avg_frame, thresh] = processFrameArray(this, rawArray, frame_range, subFrame, movieType, varargin)
            % returns a  movie using the frames specified in the input. The movie can appear different and
            % is specified with MOVIETYPE: 'orig' gives the original movie, 'diff' provides a version with a 
            % frame subtracted, by default the mean frame, 'bin' is a binary thresholded image. There is
            % an optional argument 'subFrame' specifying a frame to subtract from each frame (leave empty, [], 
            % if you don't want to specify) in order to improve the thresholding of certain objects. If giving multiple
            % thresholds, then they will be in the 4th dimension of the array.
            % RAWARRAY is in the format of dim 1,2 - height,width, 3- frame#.
            % VARARGIN{1} is the threshold level(s) for the movie, leaving
            % it out gives the default, and the returned movie is a cell
            % area of 
           
            new_mov = rawArray(:,:,frame_range);
            %detection settings
            %thresh(1) = .1; % the threshold level
            p_mouse = .0007; 
            erode_size = 3; %the size of erosion mask
            % boostContrast = 1;
            % make a movie from the average frame to subtract
            if ~isempty(subFrame)
                avg_frame = subFrame;
            else
                avg_frame = uint8(round(mean(new_mov,3)));
            end
            avg_mov = uint8(repmat(avg_frame, [1 1 size(new_mov,3)]));
            %diff_mov = imabsdiff(new_mov, avg_mov); %this should give a nice moving blob.
            diff_mov = new_mov - avg_mov; %this should give a nice moving blob.
            if(this.boostContrast)
                diff_mov = increaseMovContrast(diff_mov);
            end
            % Set the threshold for making a binary image
            thresh = this.default_thresh; %.08-.12 have worked well after image normalization
            if ~isempty(varargin) && ~isempty(varargin{1}) % I'm not exactly sure MATLAB is making empty cells
                thresh = varargin{1}; 
                if thresh == 0
                    % The way we are determining the threshold value is to take the brightest p_mouse
                    % proportion of pixels
                    sorted = sort(diff_mov(:), 'descend');
                    thresh_i = round(p_mouse*length(sorted));
                    thresh = sorted(thresh_i);
                    %thresh = .6 * graythresh(diff_mov(:)); % A way of doing auto thresholding
                    this.used_thresh = thresh;
                end
            end
            if strcmp(movieType, 'bin')
                ret_mov = [];
                for kk = 1:length(thresh)
                    bin_mov = diff_mov > thresh(kk);
                    %bin_mov = false(size(diff_mov));
                    nFrames = size(bin_mov,3);
                    for ii = 1:nFrames
                        %this is to get rid of the jagged edges due to something regarding movie compression, but also
                        %removes one pixel around the edge of the contiguous sections of blob
                        bw = bin_mov(:,:,ii);
                        %bw = im2bw(diff_mov(:,:,ii),thresh(kk));
                        bw = imerode(bw, strel('square',erode_size)); 
                        bin_mov(:,:,ii) = bw;
                    end
                    ret_mov = cat(4, ret_mov, bin_mov);
                end
            elseif strcmp(movieType, 'diff')
                ret_mov = diff_mov;
            else % show the original movie
                if(this.boostContrast)
                    new_mov = increaseMovContrast(new_mov);
                end
                ret_mov = new_mov;
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