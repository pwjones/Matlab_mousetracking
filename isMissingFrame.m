function missing = isMissingFrame(fcLum, fcPeriod, crit)
% function missing = isMissingFrame(fcLum, fcPeriod)
% 
% Tells you whether or not there is a missing frame.  Still kinda dumb
% since it doesn't tell you which one. Detects by just putting the periodic
% signal into a matrix and checking if the inter-period deviation is over a certain
% criteron amount

sig = normalizeSignal(fcLum); %gets the normalized version to work with
nPer = floor(length(fcLum)/fcPeriod);
fc2 = reshape(sig(1:(nPer*fcPeriod)), fcPeriod, []);
between = diff(fc2, 1, 2);
between_std = std(between(:));
if between_std > crit
    missing = 1;
    figure;
    plot(fc2);
    xlabel('Frame in cycle');
    ylabel('Normalized Luminance');
else
    missing = 0;
end
