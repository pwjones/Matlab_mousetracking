function fixTrackingInFolder(foldername, tracker)


vids = processVideoFolder(foldername, tracker);

s = matlabpool('size');

if s==0 %check if parallel toolbox is running.  If not, just do regular for loop
    for ii = 1:length(vids)
        vid = vids(ii);
        vid.p_mouse = 7e-4;
        vid.clearCalcData();
        vid.mousePosition([]);
        vid.save();
        vids(ii) = vid;
    end
else
    parfor ii = 1:length(vids)
        disp(['In parfor loop: ' num2str(ii)]);
        vid = vids(ii);
        vid.p_mouse = 7e-4;
        vid.clearCalcData();
        vid.mousePosition([]);
        vid.save();
        vids(ii) = vid;
    end
end