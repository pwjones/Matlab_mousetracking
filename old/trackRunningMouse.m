% script to track a mouse when running
crop = [40 40 40 40];
crop = [0 0 0 0];
mov_time_range = [16 56];
mov_time_range = [15 28];
mov_time_range = [29 110];
mov_time_range = [19.5 46];
mov_time_range = [0 29];
mt = MouseTracker('/Users/pwjones/Movies/', [], mov_time_range, crop);
track_time_range = [];
track_time_range = [11 30];
mt.mousePosition(track_time_range);

%% 
vel_range = [100 1000]; %velocity threshold for running
runFrames = find(mt.vel > vel_range(1) & mt.vel < vel_range(2));
mt.plotVelocity(runFrames);
mt.plotOrientation(runFrames);

%%
