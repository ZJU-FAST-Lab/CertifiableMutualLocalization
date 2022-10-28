
numTimes = 40; 
numRobots = 2;
numExperiments = 1;
param.noise = 0.0;
param.if_bias = true;
param.if_scale = true;
param.only_yaw = false;

lm_param.opt_using_gt = false;
lm_param.mode = "SB";
lm_param.resolution = 1;
tro_param.opt_using_gt = false;
tro_param.mode = "SB";
tro_param.resolution = 1;

for i = 1:numExperiments
    [local_trajs, bearings_data, data_gt, ~] = get_multiple_observation(numRobots, numTimes, param);
    disp("Start solving")
    %% Semidefinite Relaxation
    [output] = semidefinite_relaxation(local_trajs, bearings_data, data_gt, param);
    %% LM
    [output] = benchmark_lm(local_trajs, bearings_data, data_gt, lm_param);
    %% AM
    [output] = benchmark_tro(local_trajs, bearings_data, data_gt, tro_param);
end


%% Initial Value Sensitivity Experiment
lm_param.resolution = 4;
tro_param.resolution = 4;
[local_trajs, bearings_data, data_gt, ~] = get_multiple_observation(numRobots, numTimes, param);
[output] = benchmark_lm(local_trajs, bearings_data, data_gt, lm_param);
figure(1)
imagesc(output.cost_table)
[output] = benchmark_tro(local_trajs, bearings_data, data_gt, tro_param);
figure(2)
imagesc(output.cost_table)


