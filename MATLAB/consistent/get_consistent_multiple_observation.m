function [local_trajs, bearings_data, data_gt, world_trajs] = get_consistent_multiple_observation(numRobots, numTimes, param)
%%
% local_trajs [1, numRobots]
% bearings_data [numRobots, numRobots]
% data_gt.T [numRobots, numRobots]
% world_trajs [1, numRobots]
%%

%% Generate world_trajs
world_trajs = repmat(Trajectory(),1,numRobots);
for i = 1:numRobots
    traj = Trajectory(numTimes);
    if(param.line)
        param.x = i;
        traj.line(param);
    else
        traj.random(param);
    end
    world_trajs(i) = traj;
end

%% Generate bearings_data and T_gt
bearings_data = repmat(Bearings(),numRobots,numRobots);
T_gt = repmat(Pose(),numRobots,numRobots);
for i = 1:numRobots
    ObserverInitPose = world_trajs(i).T(1);
    ObserverInitPoseInv = Pose.inv(ObserverInitPose);
    for j = 1:numRobots
        if(i == j)
            continue;
        end
        b = Bearings(world_trajs(i), world_trajs(j), param.noise);
        bearings_data(i,j) = b;
        ObservedInitPose = world_trajs(j).T(1);
        RelativePose = ObserverInitPoseInv.mtimes(ObservedInitPose);
        T_gt(i,j) = RelativePose;
    end
end


%% Generate local_trajs
local_trajs = repmat(Trajectory(),1,numRobots);
for i = 1:numRobots
    local_trajs(i) = scale_traj(world_trajs(i),1);
end

data_gt.T = T_gt;
end