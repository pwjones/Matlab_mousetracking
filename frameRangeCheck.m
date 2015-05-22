function fr_out = frameRangeCheck(fr, mt)

if isempty(fr)
    fr_out = 1:mt.nFrames;
else
    fr_out = fr( fr > 0   &   fr <= mt.nFrames );
end