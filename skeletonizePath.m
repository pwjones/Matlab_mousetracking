function newPath = skeletonizePath(path, width, height)

im = false(height, width);
im(path.PixelIdxList) = 1;

skel_im = bwmorph(im, 'skel', inf);

newPath = path;
area = sum(skel_im(:));
pxi = find(skel_im == 1);

newPath.PixelIdxList = pxi;
newPath.Area = area;
[x, y] = ind2sub(size(skel_im), pxi);
newPath.PixelList = [y(:) x(:)];