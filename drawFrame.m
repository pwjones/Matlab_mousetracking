function frameImg = drawFrame(frameInfo, varargin)

% decide on the size of the returned img
if nargin > 1
    if nargin ~= 3
        disp('You must give both a width and a height to specify a size');
        frameImg = [];
        return;
    else
        width = varargin{1};
        height = varargin{2};
    end
else
    width = 1280;
    height = 1024;
end

frameImg = false(height, width); %the black binary image
nAreas = length(frameInfo.areas);
for ii=1:nAreas
    pl = frameInfo.areas(ii).PixelList;
    PxInd = height*(pl(:,1)-1) + (pl(:,2)-1) + 1;
    frameImg(PxInd) = 1;
end
%frameImg = imfill(frameImg, 'holes');

end