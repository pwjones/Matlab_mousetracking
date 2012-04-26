function com = track_mouse(makeVideo, timeRange)
% Mouse tracking function
% The idea is to read in a movie and generate the path of the mouse.
% makeVideo is a boolean value about whether to write an output video to 
% file or not


movie_folder = '/Users/pwjones/Documents/urbanlab/igor/mousetracking/'; %default movie folder
movie_fn = uigetfile([movie_folder '*.*']);
movie_fn = [movie_folder movie_fn];

% unfortunately, the MATLAB built-in movie processing functions seem to choke, so
% this is one downloaded from the matlab forums
[mov_struct, ~] = mmread(movie_fn,[],timeRange);

nFrames = length(mov_struct.frames);
% want to get the data into a 3D matrix - seems the best format
[binary_mov, avg_frame] = convertMovieToBinaryArray(mov_struct);

% track the position of the center of the white area over time
% preallocate AREAS
areas = regionprops(binary_mov(:,:,1), 'basic'); %built-in that gives stats on binary images
areas = repmat(areas, nFrames, 1); % make a struct arraay length of nFrames
for ii = 1:nFrames %have to work 1 frame at a time, unfortunately
    temp_reg = regionprops(binary_mov(:,:,ii), 'basic'); 
    if ~isempty(temp_reg)
        areas(ii) = regionprops(binary_mov(:,:,ii), 'basic'); 
    else 
        areas(ii) = areas(ii-1); %this will error on ii=1
    end
end
com = reshape([areas.Centroid], 2,[])'; % just an Nx2 array for the center of mass in each frame
diff_com = diff(com, 1, 1); %differences in COM for each x,y position
dists = sqrt(sum(diff_com.^2, 2)); %euclidian distances of each frame transition
[theta, rho] = cart2pol(diff_com(:,1), diff_com(:,2));
vel_pxsec = dists * mov_struct.rate; % Velocities in px/sec

if makeVideo
    useBinMovie = 0;
    vidWriter = VideoWriter('/Users/pwjones/Documents/urbanlab/igor/mousetracking/tracking.avi');
    open(vidWriter);
    figure;
    for ii=1:nFrames
        if useBinMovie
            imshow(binary_mov(:,:,ii)*255); hold on;
        else 
            imshow(mov_struct.frames(ii).cdata); hold on;
        end
        plot(com(ii,1), com(ii,2), 'r+', 'MarkerSize', 12, 'LineWidth',1);
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

% This attempts to plot a position/velocity plot on top of an image, but doesn't work because of 
% comflicting colormaps
% figure;
% img_ax = axes('Parent', gcf, 'Position', [.1 .1 .8, .8]);
% imshow(avg_frame); %hold on; %average image as background for a plotted path
% %colormap('bone');
% xlim([0 mov_struct.width]);
% ylim([0 mov_struct.height]);
% c_ax = axes('Parent', gcf, 'Position', [.1 .1 .8, .8], 'Visible', 'off');
% scatter(com(2:end,1), com(2:end,2), 6, -vel_pxsec(:), 'filled'); %plot the path of the tracked object (mouse)
% colormap(cm);
% xlim([0 mov_struct.width]);
% ylim([0 mov_struct.height]);
% 


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


