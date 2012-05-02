    
xedge=90;   %specifies cage edge to exclude experimenter's arm from analysis
Nmice=2;

% Preallocate movie space

    movieob=VideoReader('Video 10.avi');

    Nframes=movieob.NumberOfFrames;
    height = movieob.Height;
    width = movieob.Width;
    
    Nframes_est=300;    %specifies # of frames for background subtraction estimate
 
    
% Estimate background

    frames_dub=zeros(height,width);
    for aa=1:Nframes_est
        k=read(movieob,aa);
        frames_dub(:,:,aa)=mean(k,3);
        progress1 = (aa./Nframes_est)
    end
    
    avgF=uint8(mean(frames_dub,3));
    varF=uint8(std(frames_dub,0,3));
    clear frames_dub;

% Mouse tracking

    mouse(1:Nframes) = struct('cdata', zeros([],2),'xy_pos',zeros(Nmice,2));

    max=255;    % parameters for hyperbolic tangent function (contrast enhancement)
    offset=30;
    slope=0.2;
    thresh=300;

% specify portion of video to analyze    
start = 100;
stop = 500;
   
for aa=start:stop % specify which frames to analyze
    
    % get image frame, k
    k=read(movieob,aa);
    k=mean(k,3);
    k=uint8(k);
    
    
    % subtract background
    j=avgF-k;
    j=double(j);
    % subtract noisy background
        % hyperbolic tangent to brighten mouse
        k = max+tanh((j-offset)*slope)*max;
        
        % filter to smooth mouse
        h = fspecial('disk',3);
        k=imfilter(k,h);
        k(:,1:xedge)=0;
        
    % identify and keep mouse pixels
    j=k>thresh;
    [x,y]=find(k>thresh);
    mouse(1,aa).cdata=[x y];
    x=smooth(diff(sum(j)),5);
    y=smooth(diff(sum(j,2)),5);

    progress2 = (aa-start)./(stop-start)
end












