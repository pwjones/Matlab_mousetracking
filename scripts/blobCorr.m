nbi = mt2.noseblob;
noseblobs = mt2.areas(100,2);
noseblobs.PixelIdxList = 1;
otherblobs = mt2.areas(100,2);
otherblobs.PixelIdxList = 4;
for ii = 2:length(nbi)
    if ~isnan(nbi(ii))
        noseblobs(ii) = mt2.areas(ii,nbi(ii));
        rn = round(rand(1)*(mt2.nblobs(ii)-1))+1;
        while rn == nbi(ii)
            rn = round(rand(1)*(mt2.nblobs(ii)-1))+1;
        end
        otherblobs(ii) = mt2.areas(ii,rn);
    else
        noseblobs(ii) = noseblobs(1);
        otherblobs(ii) = otherblobs(1);
    end
end

fn = fieldnames(noseblobs);
nose_overlap = NaN*zeros(length(noseblobs)-1,1);
other_overlap = NaN*zeros(length(noseblobs)-1,1);
for ii = 2:length(noseblobs)
    nose_overlap(ii-1) = regionOverlap(noseblobs(ii-1), noseblobs(ii),fn);
    other_overlap(ii-1) = regionOverlap(noseblobs(ii), otherblobs(ii), fn);
end

[other_overlap_counts, ox] = hist(other_overlap,100);
[nose_overlap_counts, onx] = hist(nose_overlap,100);
figure; 
bar(ox,other_overlap_counts./max(other_overlap_counts), 'b'); hold on;
bar(onx, nose_overlap_counts./max(nose_overlap_counts), 'r');
xlim([0 1]);
