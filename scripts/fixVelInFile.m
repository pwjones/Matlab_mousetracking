function mt = fixVelInFile(pathstr, filename, tracker)
    basename = strtok(filename, '.');
    vidname = [basename '.avi']; %downside: skips other types of videos
    vidname = [pathstr filesep vidname];
    if exist(vidname, 'file')
        mt = tracker(vidname);
        np = mt.findNose(1:mt.nFrames);
        bv = mt.computeVelocity(1:mt.nFrames);
        nv = mt.noseVel(1:mt.nFrames);
        disp(['Total frames: ' num2str(mt.nFrames)]);
        disp(['Percent without nose: ' num2str(sum(isnan(np(:,1)))./mt.nFrames*100)]);
        disp(['Percent without nose velocity: ' num2str(sum(isnan(nv)./mt.nFrames*100))]);
        
    end