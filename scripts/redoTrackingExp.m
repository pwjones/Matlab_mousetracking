%redo tracking
i = 2:7;
parfor ii = i
   vid = exp.vids(ii);
   vid.fcArea =  [1221 959 1280 999];
   vid.clearCalcData();
   vid.mousePosition([]);
   vid.save;
end