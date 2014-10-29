function plotImageSection(im, rect, ax)
% function plotImageSection(im, rect, ax)
%
% rect is [left, top, right, bottom]

im2 = im(rect(2):rect(4), rect(1):rect(3));
axes(ax);
imshow(im2);
