% process files to look at sniffing rates and nose velocity

% Files and folders to process
clear perMouseData;
base_folder = VIDEO_ROOT;
mouse_names = {'19439', '21413', '971', '1080', '1527'};
%mouse_names = {'19439'};
folders = {'140401', '140404', '140409', '141016', '141017', '141020', '141021', '141022', '141024', '141027', '141028', '141029'};
nMice = length(mouse_names);
following_thresh = 20; %mm
mov_thresh = 50; %mm/sec

%% Load the data for each experiment
for ii = 1:nMice 
   SF = []; NV = []; NV_filt = []; NA = []; SF_ISI=[];% variables!!
   for jj = 1:length(folders)
        saved_file = [base_folder filesep folders{jj} filesep folders{jj} '_' mouse_names{ii} '.mat'];
        if exist(saved_file, 'file')
            disp(['Loading: ' saved_file]);
            load(saved_file);
            
            [SF_temp, NV_temp, NV_filt_temp, NA_temp, SF_ISI_temp] = correlateSniffFreqVelocityAlt(exp);
            SF = cat(1, SF, SF_temp(:)); 
            SF_ISI = cat(1, SF_ISI, SF_ISI_temp);
            NV = cat(1, NV, NV_temp(:));
            NV_filt = cat(1, NV_filt, NV_filt_temp(:));
            NA = cat(1, NA, NA_temp(:));
            
        end
   end
   % cell arrays to wrap them up!
   all_SF{ii} = SF;
   all_SF_ISI{ii} = SF_ISI;
   all_NV{ii} = NV;
   all_NV_filt{ii} = NV_filt;
   all_NA{ii} = NA;
end

%% Plotting of the relationships
sf_fig = figure;
%hist_ax = axes(); hold on;
corr_fig = figure;
pop_fig = figure; pop_ax = axes; hold on;
% stuff for trend line - a binned mean
bs = 25;
vel_bin = mov_thresh:bs:300;
low_vel = [50 100];
high_vel = [150 200];
low_rates = []; high_rates = [];
binned_freq = zeros(length(vel_bin)-1,nMice);
binned_std = zeros(length(vel_bin)-1,nMice);
for ii = 1:nMice
    figure(corr_fig);
    subplot(2, ceil(nMice/2), ii);
    hold on;

    % select parts of behavior where mouse is decelerating
    %sel = all_NA{ii} < 0;
    % select all of points
    sel = true(length(all_NA{ii}),1);
    vel_pts = all_NV_filt{ii}(sel);
    sf_pts = all_SF_ISI{ii}(sel);
    disp(sprintf('Mean Velocity: %.2f    Mean Sniff Rate: %.2f', nanmean(vel_pts), nanmean(sf_pts)));
    plot(vel_pts, sf_pts, 'k.', 'MarkerSize', 4);
     
    % screening for bad points
    toohi = find(vel_pts > 1000);
    if ~isempty(toohi)
        disp(sprintf('Mouse %s: High points %f%% - %f%% through the dataset', mouse_names{ii}, toohi(1)./length(vel_pts), toohi(end)./length(vel_pts)));
    end
    % Creating an average for various vels
    for jj = 1:(length(vel_bin)-1)
        bi = vel_pts >= vel_bin(jj) & vel_pts < vel_bin(jj+1);
        binned_freq(jj, ii) = nanmean(sf_pts(bi));
        binned_std(jj,ii) = nanstd(sf_pts(bi));
    end
    
    % compare in a couple of velocity ranges
    bi = vel_pts >= low_vel(1) & vel_pts <= low_vel(2);
    low_rates = cat(1, low_rates, sf_pts(bi));
    bi = vel_pts >= high_vel(1) & vel_pts <= high_vel(2);
    high_rates = cat(1, high_rates, sf_pts(bi));
    %[fv, xv] = ecdf(vel_pts);
    %plot(hist_ax, xv, fv);
    %vel_topi = find(fv >= .90, 1, 'first');
    %vel_top = xv(vel_topi);
    
    vel_x = vel_bin(1:end-1) + bs/2;
    %vel_x(vel_x > vel_top) = NaN; %cut out the ones that are too high to be sure about
    %high_sel = vel_x > vel_top;
    %binned_freq(high_sel,ii) = NaN;
    hold on; 
    plot(vel_x, binned_freq(:,ii), 'r', 'LineWidth',2);
    
    
    xlim([50 300]); 
    ylim([5 18]);
    xlabel('Vel (mm/sec)'); ylabel('Sniff freq (Hz)');
    title(mouse_names{ii});
    
    figure(sf_fig);
    subplot(2, ceil(nMice/2), ii);
    hold on;
    plot([0 20], [0 20], 'k--', 'LineWidth', 1);
    xlim([0 20]); ylim([0 20]);
    plot(all_SF_ISI{ii}, all_SF{ii}, 'b.', 'MarkerSize', 4);
    xlabel('ISI Sniff Freq');
    ylabel('Conv Sniff Freq');
    p = polyfit(all_SF_ISI{ii}, all_SF{ii}, 1);
    fitx = [0 5 10 18];
    fity = polyval(p,fitx);
    %plot(fitx, fity, 'r--');
    
    %addErrorBars(pop_ax, vel_x, binned_freq(:,ii), binned_std(:,ii), 'k');
    plot_err_poly(pop_ax, vel_x, binned_freq(:,ii), binned_std(:,ii), [0 0 0], [.5 .5 .5], 1);
end

figure(pop_fig);
plot(vel_x, binned_freq, 'k', 'LineWidth', 2);
xlim([50 300]); ylim([5 18]);

% Figure looking at the sniff rate for the two velocity ranges
figure;
dm = [low_rates, ones(length(low_rates),1)];
dm2 = [high_rates, 2*ones(length(high_rates),1)];
datam = cat(1, dm, dm2);
boxplot(datam(:,1), datam(:,2), 'plotstyle', 'traditional', 'notch', 'on');


