% process files to look at sniffing rates and nose velocity

% Files and folders to process
clear perMouseData;
base_folder = VIDEO_ROOT;
mouse_names = {'19439', '21413', '971', '1080'};
folders = {'140401', '140404', '140409', '141016', '141017', '141020', '141021', '141022'};
nMice = length(mouse_names);
following_thresh = 20; %mm
mov_thresh = 50; %mm/sec

%% Load the data for each experiment
for ii = 1:nMice 
   SF = []; NV = []; NV_filt = []; NA = [];% variables!!
   for jj = 1:length(folders)
        saved_file = [base_folder filesep folders{jj} filesep folders{jj} '_' mouse_names{ii} '.mat'];
        if exist(saved_file, 'file')
            load(saved_file);
            
            [SF_temp, NV_temp, NV_filt_temp, NA_temp] = correlateSniffFreqVelocity(exp);
            SF = cat(1, SF, SF_temp(:)); 
            NV = cat(1, NV, NV_temp(:));
            NV_filt = cat(1, NV_filt, NV_filt_temp(:));
            NA = cat(1, NA, NA_temp(:));
        end
   end
   % cell arrays to wrap them up!
   all_SF{ii} = SF;
   all_NV{ii} = NV;
   all_NV_filt{ii} = NV_filt;
   all_NA{ii} = NA;
end

%% Plotting of the relationships
figure;
% stuff for trend line - a binned mean
bs = 25;
vel_bin = mov_thresh:bs:400;
binned_freq = zeros(length(vel_bin)-1,nMice);

for ii = 1:nMice
    subplot(2, ceil(nMice/2), ii);
    hold on;

    % select parts of behavior where mouse is decelerating
    sel = all_NA{ii} < 0;
    
    plot(all_NV_filt{ii}(sel), all_SF{ii}(sel), 'k.', 'MarkerSize', 12);
     
    for jj = 1:(length(vel_bin)-1)
        bi = all_NV_filt{ii}(sel) >= vel_bin(jj) & all_NV_filt{ii}(sel) < vel_bin(jj+1);
        selSF = all_SF{ii}(sel);
        binned_freq(jj, ii) = nanmean(selSF(bi));
    end
    vel_x = vel_bin(1:end-1) + bs/2;
    hold on; 
    plot(vel_x, binned_freq(:,ii), 'r', 'LineWidth',1);
    
    xlim([50 300]); ylim([5 18]);
    xlabel('Vel (mm/sec)'); ylabel('Sniff freq (Hz)');
    title(mouse_names{ii});
end

figure;
plot(vel_x, binned_freq);
