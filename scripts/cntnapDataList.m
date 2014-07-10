% cntnapDataList.m
% Specifies the files used for the CNTNAP2 mouse dataset

mouse_names = {'421', '422', '424', '425'};
% The full list
folders = {'140618', '140619', '140620', '140623', '140625','140626', '140627'};
train_trials = {1:26, 1:26, 1:26, 1:26, 1:26, 1:26};
ctl_trials = {1:26, 1:26, 1:26, 1:26, 1:26, 1:26};
occr_trials = {[], [], [], [], [], [], [], []}; %each cell is a mouse
occl_trials = {[], [], [], [], [], [], [], []};
nMice = length(mouse_names);