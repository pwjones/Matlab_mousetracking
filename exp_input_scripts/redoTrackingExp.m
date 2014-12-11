%redo tracking
i = 1:length(exp.vids);
parfor ii = i
   vid = exp.vids(ii);
   vid.fcArea =  [1221 959 1280 999];
   vid.clearCalcData();
   vid.mousePosition([]);
   vid.save;
end

save('141028_1080.mat', 'exp');