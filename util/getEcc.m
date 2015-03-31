function ecc = getEcc(mt, frames)
% function ecc = getEcc(mt, frames)
% 
%
ecc = zeros(length(frames), mt.maxBlobs);
for ii = 1:length(frames)
    fi = frames(ii);
    areas = mt.areas(fi,:);
    ecc_temp = [areas.MajorAxisLength]./[areas.MinorAxisLength];
    ecc(ii, 1:length(ecc_temp)) = ecc_temp;
end