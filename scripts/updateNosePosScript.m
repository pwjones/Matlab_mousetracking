% complete tracking

for ii = 24:length(folders)
    vids = processVideoFolder(folders{ii},@MouseTrackerKF);
    for jj = 1:length(vids)
        for kk = 1:vids(jj).nFrames
            if ~isnan(vids(jj).noseblob(kk))
                vids(jj).nosePos(kk,:) = vids(jj).areas(kk,vids(jj).noseblob(kk)).Centroid;
            end
        end
        vids(jj).save;
    end 
end