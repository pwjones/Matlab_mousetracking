function dists = cartDist(pts)
% function dists = cartDist(pts)

dists = zeros(size(pts,1), size(pts,1), 2);
for ii = 1:2
    dists(:,:,ii) = eye(size(pts,1)) * pts(:,ii);
end