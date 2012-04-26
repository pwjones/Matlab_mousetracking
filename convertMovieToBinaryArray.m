function [binArray, avg_frame] = convertMovieToBinaryArray(mov_struct) %#codegen

% declare rgb2gray extrinsic, meaning that it stays a matlab function
coder.extrinsic('rgb2gray', 'imabsdiff', 'graythresh', 'im2bw');

mov = mov_struct.frames;
new_mov = zeros([size(mov(1).cdata,1) size(mov(1).cdata,2) length(mov)], 'uint8');
 for ii=1:length(mov)
     new_mov(:,:,ii) = rgb2gray(mov(ii).cdata);
 end
% make a movie from the average frame to subtract
avg_frame = uint8(round(mean(new_mov,3)));
avg_mov = repmat(avg_frame, [1 1 length(mov)]);
mean_sub_mov = zeros(size(new_mov)); % need to declare the array so that the MATLAB lets me index it later
mean_sub_mov = imabsdiff(new_mov, avg_mov); %this should give a nice moving blob.

binArray = zeros(size(mean_sub_mov), 'uint8');

nFrames = 1;
thresh = graythresh(mean_sub_mov);
nFrames = int16(size(binArray,3));
for ii = int16(1):nFrames
    binArray(:,:,ii) = im2bw(mean_sub_mov(:,:,ii), thresh);
end
%binArray = binArray * 255;
return;