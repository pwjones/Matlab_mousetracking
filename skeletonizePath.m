function newPath = skeletonizePath(path, width, height)

im = false(height, width);
im(path.PixelIdxList) = 1;

skel_im = bwmorph(im, 'skel', inf);

% figure; imshow(im);
% figure; imshow(skel_im);

newPath = path;
area = sum(skel_im(:));
pxi = find(skel_im == 1);

newPath.PixelIdxList = pxi;
newPath.Area = area;
x = floor(pxi/height);
y = mod(pxi,height);
newPath.PixelList = [x(:) y(:)];