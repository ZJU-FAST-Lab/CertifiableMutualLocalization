function [result] = local_optimization_lm(data)

    auxiliary_data = data.auxiliary_data;
    numRobots = data.numRobots;
    %% Run SE-Sync
    t_start = tic();

    obj3=@(x1)cost_function(x1,data);
    x1 = zeros(3*numRobots,1);
    for i = 1:numRobots
        rot = rnd.rot();
        [r,p,y] = bot_quat_to_roll_pitch_yaw(bot_matrix_to_quat(rot));
        x1((i-1)*3+1:(i-1)*3+3, 1)= [r;p;y];
    end

    [x,resnorm,residual,exitflag,output] = lsqnonlin(obj3,x1);
    R_opt = zeros(3, 3*numRobots);
    for k = 1:numRobots
        R = bot_quat_to_matrix(bot_roll_pitch_yaw_to_quat( x((k-1)*3+1:(k-1)*3+3, 1) ) );
        R_opt(:, 3*(k-1) + 1 : 3*(k-1) + 3) = R;
    end
    t_end = toc(t_start);
    g_RTR = evaluate_objective_WYJ(R_opt, auxiliary_data);
    f_X = evaluate_objective_WYJ_by_X(R_opt'*R_opt, auxiliary_data);

    result.R_opt = R_opt(:,1:3)'*R_opt;
    result.f_X = f_X;
    result.g_RTR = g_RTR;
    result.refuse_error = 0.0;
    result.rounding_error = 0.0;
    result.time = t_end;
    result.l2error = norm(result.R_opt - data.data_gt.R_gt);
    result.X = R_opt'*R_opt;
end

function [cost] = cost_function(x, data)
    auxiliary_data = data.auxiliary_data;
    numRobots = data.numRobots;
    R_opt = zeros(3, 3*numRobots);
    for k = 1:numRobots
        R = bot_quat_to_matrix(bot_roll_pitch_yaw_to_quat( x((k-1)*3+1:(k-1)*3+3, 1) ) );
        R_opt(:, 3*(k-1) + 1 : 3*(k-1) + 3) = R;
    end
    cost = evaluate_objective_WYJ(R_opt, auxiliary_data);
end