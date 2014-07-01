% Concentration variation data

mouse_names = {'19439-1', '19439-2', '19439-3', '21413-1', '21413-2', '21413-3'};
% The full list
folders = {'140610', '140611', '140612', '140613', '140616','140617', '140618'};
train_trials = {1:21, 1:21, 1:21, 1:21, 1:21,1:21};
ctl_trials = {1:21, 1:21, 1:21, 1:21, 1:21, 1:21};
occr_trials = {[], [], [], [], [], [], [], []}; %each cell is a mouse
occl_trials = {[], [], [], [], [], [], [], []};
nMice = length(mouse_names);