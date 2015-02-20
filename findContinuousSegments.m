function ranges = findContinuousSegments(inds)

jumps = diff(inds); %identify contiguous and non are 
skips = find(jumps > 1); % noncontiguous
if ~isempty(inds)
    seq_end = [inds(skips)' inds(end)]'; %these are the indices of segments of following
    seq_start = [inds(1) inds(skips+1)']';
    if (isempty(seq_end)) seq_end = inds(end); end
    ranges = [seq_start(:) seq_end(:)];
else 
    ranges = [];
end