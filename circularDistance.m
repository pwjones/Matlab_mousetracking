function cdist = circularDistance(a1, a2)
% finds the phase difference of angles a1 and a2 (expressed in radians).  
% a1 is the reference, so the difference is a1-a2
if isempty(a1) || isempty(a2)
    cdist = [];
else
    cdist = mod(a1-a2, 2*pi);
    ap = cdist>pi;
    cdist(ap) = 2*pi-cdist(ap);
end
