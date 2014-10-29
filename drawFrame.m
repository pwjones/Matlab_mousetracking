function frameImg = drawFrame(frameInfo)

width = 1280;
height = 1024;

frameImg = false(height, width); %the black binary image
nAreas = length(frameInfo.areas);
for ii=1:nAreas
    pl = frameInfo.areas(ii).PixelList;
    PxInd = height*(pl(:,1)-1) + (pl(:,2)-1) + 1;
    frameImg(PxInd) = 1;
end

end