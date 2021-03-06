function [paraRel, absAng] = convertMouseHeadingAngles(orthoRel, ortho, pos)
% function para_angles = convertMouseHeadingAngles(ortho_angles, pos)
% Output: paraRel - The angle relative to the trail angle - 0 is parallel
%         AbsAngle - The absolute angle relative to the horizontal
%
% The purpose of this is to convert from Orthogonal angles to parallel 
% angles.
absAng = orthoRel + ortho;
isright = pos > 0; %orthogonal distance from trail. Right is positive by convention
para = ortho; 
para(isright) = ortho(isright) + pi/2; 
para(~isright) = ortho(~isright) + 3*pi/2;
paraRel = circ_dist(absAng, para); 