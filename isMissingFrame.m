function missing = isMissingFrame(fcLum, fcPeriod, crit)
% function missing = isMissingFrame(fcLum, fcPeriod)
% 
% Tells you whether or not there is a missing frame.  Still kinda dumb
% since it doesn't tell you which one. Detects by just putting the periodic
% signal into a matrix and checking if the inter-period deviation is over a certain
% criteron amount
pb = 1;

sig = normalizeSignal(fcLum); %gets the normalized version to work with
nPer = floor(length(fcLum)/fcPeriod);
fc2 = reshape(sig(1:(nPer*fcPeriod)), fcPeriod, []);
between = diff(fc2, 1, 2);
between_std = std(between(:));
mean_lum = mean(fc2, 2);

if(pb)
    figure; hold on;
    plot(mean_lum, 'LineWidth',4);
    plot(fc2);
    xlabel('Frame in cycle');
    ylabel('Normalized Luminance');
end
missing = between_std > crit;

