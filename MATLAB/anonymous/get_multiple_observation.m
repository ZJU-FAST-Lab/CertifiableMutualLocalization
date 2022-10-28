function [local_trajs, bearings_data, data_gt, world_trajs] = get_multiple_observation(numRobots, numTimes, param)
%%
% local_trajs [1, numRobots]
% bearings_data [1, numRobots-1]
% world_trajs [1, numRobots-1]
%%

%% Generate world_trajs
world_trajs = repmat(Trajectory(),1,numRobots);
for r = 1:numRobots
    traj = Trajectory(numTimes);
    traj.random(param);
    world_trajs(r) = traj;
end

%% Generate bearings_data and T_gt
bearings_data = repmat(Bearings(),1,numRobots-1);
T_gt = repmat(Pose(),1,numRobots-1);
ObserverInitPose = world_trajs(1).T(1);
ObserverInitPoseInv = Pose.inv(ObserverInitPose);
for r = 2:numRobots
    bearings = Bearings(world_trajs(1), world_trajs(r), param.noise);
    bearings_data(r-1) = bearings;
    
    ObservedInitPose = world_trajs(r).T(1);
    RelativePose = ObserverInitPoseInv.mtimes(ObservedInitPose);
    T_gt(r-1) = RelativePose;
end

%% Generate local_trajs
local_trajs = repmat(Trajectory(),1,numRobots);
local_trajs(1) = scale_traj(world_trajs(1),1);
scale_gt = zeros(1,numRobots-1);
for r = 2:numRobots
    if(param.if_scale)
        scale = rand*10;
    else
        scale = 1.0;
    end
    local_trajs(r) = scale_traj(world_trajs(r),scale);
    scale_gt(r-1) = 1.0 / scale;
end

order_GT = randperm(numRobots-1);
bearings_data = bearings_data(:,order_GT);

%% Pack GroundTruth
data_gt.order = order_GT;
data_gt.T = T_gt;
data_gt.s = scale_gt;
end