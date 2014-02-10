mt.blob_num = 1;
for ii=1:mt.nFrames
    mt.assignBlobIDs(ii);
end

overlaps = nan*zeros(mt.nFrames, 1);
for ii=2:mt.nFrames
    if ~isnan(mt.noseblob(ii)) && ~isnan(mt.noseblob(ii-1))
        overlaps(ii) = regionOverlap(mt.areas(ii, mt.noseblob(ii)), mt.areas(ii-1, mt.noseblob(ii-1)), mt.regprops);
    end
end

    
for ii=1:mt.nFrames
    if ~isnan(mt.noseblob(ii))
        mt.nosePos(ii,:) = mt.areas(ii,mt.noseblob(ii)).Centroid;
    end
end
