function com = track_mouse(makeVideo, timeRange)
% Mouse tracking function
% The idea is to read in a movie and generate the path of the mouse.
% makeVideo is a boolean value about whether to write an output video to 
% file or not

MAX_FRAMES = 400; % the maximum number of frames processable at one time by matlab due to memory limitations

movie_folder = '/Users/pwjones/Documents/urbanlab/data/mouse_training/'; %default movie folder
movie_fn = uigetfile([movie_folder '*.*']);
if movie_fn == 0 % the user has canceled the file selection
    com = [];
    return;
end
movie_fn = [movie_folder movie_fn];

% unfortunately, the MATLAB built-in movie processing functions seem to choke, so
% this is one downloaded from the matlab forums
disp('Reading movie');
[mov_struct, ~] = mmread(movie_fn,[],timeRange);

disp('Converting to movie to grayscale');
mov_struct = convertToGray(mov_struct); %converting to gray immediately is probably good.

% want to get the data into a 3D matrix - seems the best format.  
% However, there is some complication with large movies fitting into memory.
% Need to divide the movie into segments, and process each sequentially.
% Can not parallelize on a single computer because this is a memory issue.
nFrames = length(mov_struct.frames);
nSegs = ceil(nFrames/MAX_FRAMES);
binary_mov = zeros([mov_struct.height, mov_struct.width, nFrames], 'uint8');
%computes the average frame over the whole movie using a subset of frames
avg_frame = averageFrame(mov_struct, 1:10:nFrames); 
%avg_frame = zeros([mov_struct.height, mov_struct.width, nSegs], 'uint8');
for ii=1:nSegs
    disp(['Processing segment ' num2str(ii)]);
    first = (ii-1)*MAX_FRAMES + 1; %first frame of segment
    if (ii == nSegs) % the last frame
        last = nFrames;
    else
        last = (ii)*MAX_FRAMES;
    end
        [binary_mov(:,:,first:last), ~] = convertMovieToBinaryArray(mov_struct, first:last, avg_frame);
end

%avg_frame = uint8(mean(avg_frame,3));

% track the position of the center of the white area over time
% preallocate AREAS
% List of attributes to calculate about the detected areas
props = {'BoundingBox', 'Centroid', 'Perimeter', 'Area', 'Orientation', 'MajorAxisLength', 'MinorAxisLength'};
areas = regionprops(binary_mov(:,:,1), props); %built-in that gives stats on binary images
areas = repmat(areas, nFrames, 1); % make a struct arraay length of nFrames
for ii = 1:nFrames %have to work 1 frame at a time, unfortunately
    temp_reg = regionprops(binary_mov(:,:,ii), props); 
    if ~isempty(temp_reg)
        areas(ii) = regionprops(binary_mov(:,:,ii), props); 
    else 
        areas(ii) = areas(ii-1); %this will error on ii=1
    end
end
% Centers of mass (COM)
com = reshape([areas.Centroid], 2,[])'; % just an Nx2 array for the center of mass in each frame
% Compute the velocities
diff_com = diff(com, 1, 1); %differences in COM for each x,y position
dists = sqrt(sum(diff_com.^2, 2)); %euclidian distances of each frame transition
[theta, rho] = cart2pol(diff_com(:,1), diff_com(:,2));
vel_pxsec = dists * mov_struct.rate; % Velocities in px/sec
% Bounding boxes
bounds = reshape([areas.BoundingBox], 4, [])'; %reshape to Nx4 array
% Shape based orientations
theta = [areas.Orientation];
% ellipse measurements
maj_ax_l = [areas.MajorAxisLength];
min_ax_l = [areas.MinorAxisLength];

if makeVideo
    useBinMovie = 1;
    vidWriter = VideoWriter([movie_folder '/tracking.avi']);
    open(vidWriter);
    figure;
    for ii=1:nFrames
        if useBinMovie
            imshow(binary_mov(:,:,ii)*255); hold on;
        else 
            imshow(mov_struct.frames(ii).cdata); hold on;
        end
        plot(com(ii,1), com(ii,2), 'r+', 'MarkerSize', 12, 'LineWidth',1);
        %rectangle('Position',bounds(ii,:), 'Curvature', [1 1], 'EdgeColor', 'r', 'FaceColor', 'none');
        ellipse(areas(ii).MajorAxisLength,areas(ii).MinorAxisLength, theta(ii)./180*pi, com(ii,1),com(ii,2),'r');
        
        currFrame = getframe;
        writeVideo(vidWriter,currFrame);
        hold off;
    end
    close(vidWriter);
end

% plotting the path over the average frame
figure;
cm = colormap(jet(nFrames));
sorted_vel = sort(vel_pxsec); %sorted vector of speeds for colormapping
imshow(avg_frame); hold on;
xlim([0 mov_struct.width]);
ylim([0 mov_struct.height]);
% unfortunately, I haven't been able to figure out an easier way to do colormapping of points plotted over 
% a B&W image, without the points and the image sharing the same colormap.  So, here I'm just plotting every
% point separately, with a color indicative of the velocity of the animal.  
for ii=2:nFrames
    ci = find(sorted_vel == vel_pxsec(ii-1), 1, 'first');
    line('Xdata', com(ii,1), 'Ydata', com(ii,2), 'Marker', '.', 'MarkerSize', 10, 'Color', cm(ci,:));
end


% plotting the position/direction (based on movement)
figure;
cm = colormap(hsv(nFrames));
sorted_dir = sort(theta); %sorted vector of speeds for colormapping
imshow(avg_frame); hold on;
xlim([0 mov_struct.width]);
ylim([0 mov_struct.height]);
% unfortunately, I haven't been able to figure out an easier way to do colormapping of points plotted over 
% a B&W image, without the points and the image sharing the same colormap.  So, here I'm just plotting every
% point separately, with a color indicative of the velocity of the animal.  
for ii=2:nFrames
    ci = find(sorted_dir == theta(ii-1), 1, 'first');
    line('Xdata', com(ii,1), 'Ydata', com(ii,2), 'Marker', '.', 'MarkerSize', 10, 'Color', cm(ci,:));
end

% -----------------------------------------------
function mov_struct = convertToGray(mov_struct)
% function that converts the movie, frame by frame into grayscale - if it's not already.
%

for ii=1:length(mov_struct)
    mov_struct.frames(ii).cdata = rgb2gray(mov_struct.frames(ii).cdata);
end

% ---------------------------------------------------------
function avg_frame = averageFrame(mov_struct, frame_range) 
% function to that takes the specific movie structure and computes the average frame of the specified frames

mov = mov_struct.frames(frame_range);
new_mov = zeros([size(mov(1).cdata,1) size(mov(1).cdata,2) length(mov)], 'uint8');
for ii=1:length(mov)
    new_mov(:,:,ii) = squeeze(mov(ii).cdata(:,:,1));
end
avg_frame = uint8(round(mean(new_mov,3)));

% --------------------------------------------------------------------------------------------
function [bin_mov, avg_frame] = convertMovieToBinaryArray(mov_struct, frame_range, subFrame)
% returns a binary thresholded movie using the frames specified in the input. Also there is 
% an optional argument 'subFrame' specifying a frame to subtract from each frame in order to 
% improve the thresholding of certain objects.

mov = mov_struct.frames(frame_range);
new_mov = zeros([size(mov(1).cdata,1) size(mov(1).cdata,2) length(mov)], 'uint8');
%new_mov = [mov.cdata];
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

bin_mov = zeros(size(new_mov), 'uint8');
nFrames = size(bin_mov,3);
thresh = graythresh(new_mov);
for ii = 1:nFrames
    %thresh = graythresh(new_mov(:,:,ii));
    bin_mov(:,:,ii) = im2bw(new_mov(:,:,ii), thresh);
end

% --------------------------------------------------------
function com = findCOM(frame)
% This is a function that computes the center of mass of the white pixels in a binary image
% Unfortunately, currently RETURNS INCORRECT X-VALUES. Need to FIX.  
thresh = 100;
f_size = size(frame); 
white_px = find(frame > thresh);
y = mod(white_px, f_size(1));
x = white_px / f_size(2);
%x2 = ((white_px-1) ./ f_size(2)) + 1;
com = [mean(x) mean(y)];


