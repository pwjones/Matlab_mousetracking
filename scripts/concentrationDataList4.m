% These are two experienced mice tested with Beta-ionone under randomly interleaved trial conditions,
% where care was taken to normalize the difficulty of the trails between concentrations in terms of
% general curviness and number of crossings and near-crossings.  
% Data was collected starting on 7/18/2014 till 
% Data
mouse_names = {'19439-00005', '19439-0001', '19439-001', '19439-01', '21413-00005', '21413-0001', '21413-001', '21413-01'};
folders = {'140718', '140721','140722', '140723', '140724', '140725'};
train_trials = {1:6, 1:10, 1:10, 1:10, 1:5, 1:10, 1:10, 1:10};
ctl_trials = {1:6, 1:10, 1:10, 1:10, 1:5, 1:10, 1:10, 1:10};
occr_trials = {[], [], [], [], [], [], [], []}; %each cell is a mouse
occl_trials = {[], [], [], [], [], [], [], []};
nMice = 2;
nConc = 4;