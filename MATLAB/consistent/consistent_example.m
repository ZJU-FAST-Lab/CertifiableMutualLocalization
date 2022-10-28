%% Reset environment
clear all;
close all;
clc;

numRobots = 2;
numExperiments = 1;
numTimes = 40;
param.noise = 0.0;
param.if_bias = false;
param.if_scale = false;
param.only_yaw = false;
param.line = false;
param.mode = "CrossProduct"; %"Marginalization"
for indexExperiment = 1:numExperiments
    [local_trajs, bearings_data, data_gt, world_trajs] = get_consistent_multiple_observation(numRobots, numTimes, param);
    [auxiliary_data, data_gt] = generate_input(local_trajs, bearings_data, world_trajs, param);
    
    input.numRobots = numRobots;
    input.auxiliary_data = auxiliary_data;
    input.data_gt = data_gt;
    input.tolerance = 1e-1;
    
    %% Iteration
    iteration_result = convex_iteration(input);
    
    %% SDP
    sdp_result = pure_sdp(input);

    %% Local Optimization
    local_result = local_optimization(input);

    %% Local Optimization_lm
    lm_result = local_optimization_lm(input);
    
    %% GroundTruth
    f_X_gt = evaluate_objective_WYJ_by_X(data_gt.Z_gt, auxiliary_data);
    g_RTR_gt = evaluate_objective_WYJ(data_gt.R_gt, auxiliary_data);
    [V,D] = eig(data_gt.Z_gt); U = V(:,1:end-3); C = U * U';
    f_X_CX_gt = f_X_gt + trace(C*data_gt.Z_gt);
    
    disp(sprintf('g(RTR): SDP = %g,  Iteration = %g, GT = %g, local = %g\n' , sdp_result.g_RTR, iteration_result.g_RTR, g_RTR_gt, local_result.g_RTR));
    disp(sprintf('f(X): SDP = %g,  Iteration = %g, GT = %g, local = %g\n' , sdp_result.f_X, iteration_result.f_X, f_X_gt, local_result.f_X));
    disp(sprintf('f(X) + trace(CX): Iteration = %g,  GT = %g\n', iteration_result.f_X_CX, f_X_CX_gt));
    disp(sprintf('residual_norm: SDP = %g,  Iteration = %g, local = %g\n', sdp_result.l2error, iteration_result.l2error,local_result.l2error));
    disp(sprintf('time: SDP = %g,  Iteration = %g, local = %g\n', sdp_result.time, iteration_result.time, local_result.time));

   

    R_opt = iteration_result.R_opt;
    %% Get Distance
    distance_b = auxiliary_data.distance_b_static;
    distance_b_wait = auxiliary_data.distance_b_wait;
    distance_A = auxiliary_data.distance_A;
    dim_distance_variable = auxiliary_data.dim_distance_variable;
    for i = 1:dim_distance_variable
        distance_b(i,1) = distance_b(i,1) + trace(R_opt' * R_opt * distance_b_wait{i,1});
    end
    [L,U] = lu(distance_A);
    distance_A_inv_lu = inv(U) * inv(L);
    
    distance_opt = - distance_A_inv_lu * distance_b;
    if(sum(distance_opt) < 0)
        distance_opt = - distance_opt;
    end
    l2error_distance = norm(distance_opt  - data_gt.distance_gt)


    %% Get Translation
    T_var = zeros(3*numRobots+1, 3*numRobots+1);
     for i = 1:numRobots
        R_i = R_opt(:,3*(i-1)+1:3*i);
        for j = 1:numRobots
            R_j = R_opt(:,3*(j-1)+1:3*j);
            if(i==j)
                continue
            end
            for t = 1:numTimes
                R_i_t = local_trajs(i).T(t).R; t_i_t = local_trajs(i).T(t).t;
                R_j_t = local_trajs(j).T(t).R; t_j_t = local_trajs(j).T(t).t;

                b_i_t = bearings_data(i,j).bearing(t).x; d_i_t = bearings_data(i,j).depth(t);

                if(j > i)
                    virtual_j = j-1;
                else
                    virtual_j = j;
                end
                index_distance_i_t = (i-1)*(numRobots-1)*numTimes + (virtual_j-1)*numTimes + t;
                distance = distance_opt(index_distance_i_t, 1);
                C_var = zeros(3, 3*numRobots+1);
                C_var(:,3*(i-1)+1 : 3*(i-1)+3) = eye(3);
                C_var(:,3*(j-1)+1 : 3*(j-1)+3) = -eye(3);
                C_var(:,end) = R_i* (R_i_t * b_i_t * d_i_t + t_i_t) - R_j * t_j_t;
                T_var = T_var + C_var' * C_var;
            end
        end
     end

    A = T_var(1:end-1, 1:end-1);
    b = T_var(1:end-1, end);
    vect_opt = -pinv(A) * b;
    first_t_opt = vect_opt(1:3,1);
    t_opt = zeros(3,numRobots);
    for i = 1:numRobots
        t_opt(:,i) = vect_opt(3*(i-1)+1:3*(i-1)+3, 1) - first_t_opt;
    end
    
    t_gt = data_gt.t_gt;
    first_t_gt = t_gt(:,1);
    for i = 1:numRobots
        t_gt(:,i) = t_gt(:,i) - first_t_gt;
    end
    l2error_translation = norm(t_opt  - t_gt);
end






