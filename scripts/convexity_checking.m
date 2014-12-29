for ii = 1:length(a)
    elongation(ii) = a(ii).MajorAxisLength / a(ii).MinorAxisLength;
    ratio(ii) = a(ii).Perimeter./sqrt(a(ii).Area);
    pToBound(ii) = a(ii).Perimeter ./ (2*a(ii).BoundingBox(3) + 2*a(ii).BoundingBox(4));
    areaToEArea(ii) = 1-(a(ii).Area / (a(ii).MajorAxisLength/2 * a(ii).MinorAxisLength/2 * pi));
end
elongation
ratio
pToBound
areaToEArea
