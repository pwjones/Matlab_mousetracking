% Concentration variation data

mouse_names = {'19439-1', '19439-2', '19439-3', '21413-1', '21413-2', '21413-3'};
% The full list
folders = {'140521', '140522', '140603', '140604', '140605', '140609'};
train_trials = {1:17, 1:17, 1:17, 1:17, 1:17,1:16};
ctl_trials = {1:17, 1:17, 1:17, 1:17, 1:17, 1:16};
occr_trials = {[], [], [], [], [], [], [], []}; %each cell is a mouse
occl_trials = {[], [], [], [], [], [], [], []};
nMice = length(mouse_names);