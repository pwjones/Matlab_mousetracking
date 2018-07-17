%Generate fake following data
% all measurements in mm, mm/msec, radians
% strategy would be some random noise in the direction w/ some purposeful
% turning at random times
% Positions are in distance from the simulated trail
simFollowing = [];
for jj=1:100
    ntp = 100;
    dt = 20; %ms
    %vel = 60 + normrnd(40, 20, ntp, 1);
    vel = 100 / 1000; % mm/msec
    len = 100;
    %mu = normrnd(0, pi/4, ntp, 1);
    %mu = random('unif', -pi/2, pi/2, [ntp, 1]);
    %mu = random('norm', 0, pi/8, [ntp, 1]);
    sign = 2*(rand(ntp,1) > .5) - 1;
    mu = random('burr', 2.36, 0.974, 5.251, [ntp, 1]) .* sign; %fit to turn data
    theta = cumsum(mu);
    y = zeros(ntp, 1);
    x = zeros(ntp, 1);
    for ii=2:ntp
        y(ii) =y(ii-1) + (vel*dt)*sin(theta(ii));
        x(ii) = x(ii-1) + (vel*dt)*cos(theta(ii));
    end
    trail = [0:(len-1)]';
    trail = [trail zeros(len,1)];
    dists = orthogonalDistance([x,y], trail, theta);
    outside = find(dists > 20 | dists < -20, 1, 'first');
    simFollowing = cat(1, simFollowing, dists(1:outside));
end