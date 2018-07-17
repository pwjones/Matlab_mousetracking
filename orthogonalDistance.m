function [dists, vects] = orthogonalDistance(nosePos, trailPos, headingTheta)
%
% nosePos - the first vector of positions
% trailPos

%nosePos = this.nosePos(frames,:);
%trailPos = this.pathVertices(pathNum).PixelList;
distm = ipdm(single(nosePos), single(trailPos));

% Must figure out if the animals' nose is over the trail or
% not. Do this by getting the 4
% closest trail points to each nose position and see if they encircle it.
noseDist = NaN*zeros(size(nosePos,1),6); 
closestTrailP = ones(size(nosePos,1),2,6);
for ii = 1:6
    [noseDist(:,ii), mini] = nanmin(distm, [], 2); %get the minimum value
    li = sub2ind(size(distm), (1:size(nosePos,1))', mini);
    distm(li) = NaN; %set the mins to NaN to get the next closest on next iteration
    closestTrailP(:,:,ii) = trailPos(mini, :);
end
over = isContained(nosePos, closestTrailP);

noseDist = noseDist(:,1);
vects = closestTrailP(:,:,1) - nosePos;
noseDist(over) = 0; %set those points to zero because the nose IS over the trail itself
% Now we use the heading of the vector to give the direction from the trail
% a sign
headingTheta = mod(headingTheta, 2*pi);
orthoTheta = mod(cart2pol(vects(:,1), vects(:,2)),2*pi); %the orientation of the vector to the closest trail point
[rotations, ~] = rotationDirection(headingTheta, orthoTheta);
dists = noseDist;
dists = dists .* rotations; %gives it a sign