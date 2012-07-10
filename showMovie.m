function showMovie(seq)

global exitmovie;

exitmovie = 0;
fh = figure;

set(fh, 'WindowKeyPressFcn', @exitMovieLoop);
nframes = size(seq,3);
for ii = 1:nframes
    imshow(seq(:,:,ii));
    pause(1/60);
    if exitmovie
        break;
    end
end


function exitMovieLoop(this, src, event)
% Function to set a flag internally to exit a movie that is being displayed.
% It is not used for any other purpose

global exitmovie;

if strcmp(event.Key, 'q') || strcmp(event.Key, 'escape')
    exitmovie = 1;
end