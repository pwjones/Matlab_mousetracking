function [p,h] = testCellArrayMedians(c1, c2)

a1 = []; %gonna just make a long array
a2 = [];

for ii = 1:length(c1)
    a1 = cat(1, a1, c1{ii}(:));
end

for ii = 1:length(c2)
    a2 = cat(1, a2, c2{ii}(:));
end

[p,h] = ranksum(a1, a2);
