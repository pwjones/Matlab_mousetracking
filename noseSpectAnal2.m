function [f, S, tah, fah] = noseSpectAnal2(traj, dt, varargin)

s = size(traj);
dbg = 0;
%ntrial = s(2);
if nargin > 2
    range = varargin{1};
else 
    range = [0 1/dt/2];
end
param.tapers = [2 3]; %[Time-Bandwidth Product (temporal x spectrum variance of pulses), N tapers (2TWP-1)]
param.Fs = 1/dt; %sampling frequency
param.pad = 0; % padding option, -1=none, 0=next pow2, 1= next,next pow 2, ....
range_high = min(param.Fs/2, range(2));
param.fpass = [range(1) range_high]; % Frequency range returned
param.err = 0; %provide error estimates?
param.trialave = 1; %average trials?

[S,f] = mtspectrumc(traj, param);
%[S,f] = mtspectrumsegc(traj, 5, param, 1);

tah = []; fah = [];
if dbg
    
    figure;
    t = (0:(s(1)-1)) * dt;
    tah = subplot(2,1, 1);
    ntrace = min(10, s(2));
    plot(t, traj(:,1:ntrace));
    title('Nose position (1000 frames)')
    xlabel('time (seconds)');
    ylabel('Distance (mm)');
    
    fah = subplot(2,1,2);
    plot(f,S);
    xlabel('Frequency (Hz)');
    ylabel('Power');
end
