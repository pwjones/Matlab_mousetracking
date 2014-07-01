figure; hold on;
bodyVel = []; noseVel = [];
for jj = 1:length(perMouseData)
    for kk = 1:length(perMouseData(jj).body_vel)
        for ll= 1:length(perMouseData(jj).body_vel{kk})
            temp = perMouseData(jj).body_vel{kk}{ll};
            bodyVel = cat(1, bodyVel, temp(:));
            temp = perMouseData(jj).nose_vel{kk}{ll};
            noseVel = cat(1, noseVel, temp(:));
        end
    end
end