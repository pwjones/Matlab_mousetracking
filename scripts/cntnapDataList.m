% cntnapDataList.m
% Specifies the files used for the CNTNAP2 mouse dataset

mouse_names = {'421', '422', '424', '425'};
genotype = [1 2 2 3]; %1 - WT, 2-Het, 3-Homo
gcolor = {[0 0 0], [.1 .2 .9], [.6 .7 1]};
% The full list
folders = {'140618', '140619', '140620', '140623', '140625','140626', '140627', '140630', ...
           '140701', '140702', '140703', '140707', '140708'};
train_trials = {1:53, 1:53, 1:54, 1:53, 1:53, 1:53};
ctl_trials = {1:53, 1:53, 1:54, 1:53, 1:53, 1:53};
occr_trials = {[], [], [], [], [], [], [], []}; %each cell is a mouse
occl_trials = {[], [], [], [], [], [], [], []};
nMice = length(mouse_names);