function  [output]  = benchmark_lm(local_trajs, bearings_data, data_gt, lm_param)
%%
% local_trajs [1, numRobots]
% bearings_data [1, numRobots-1]
%%

T_gt = data_gt.T(1);
R_gt = T_gt.R;
resolution = lm_param.resolution;
mode = lm_param.mode;
opt_using_gt = lm_param.opt_using_gt;
numTimes = size(local_trajs(1).poses,2);

output.cost_table = zeros(resolution, resolution);
output.times = zeros(resolution, resolution);
output.iteration = zeros(resolution, resolution);
output.l2error_R = zeros(resolution, resolution);
output.l2error_t = zeros(resolution, resolution);

obj3=@(x_guess)costMarg(x_guess,local_trajs,bearings_data,mode);

for i = 1:resolution
    roll = i * 2 * pi / resolution - pi;
    for j = 1:resolution
        [i,j]
        pitch = j * 2 * pi/ resolution- pi;
        rpy = [roll; pitch; 0];
        perturbance = bot_quat_to_matrix(bot_roll_pitch_yaw_to_quat(rpy));
        
        x0 = perturbance*R_gt;
        if(resolution == 1)
            x0 = rnd.rot();
        end
        if(opt_using_gt)
            x0 = R_gt;
        end

        [r,p,y] = bot_quat_to_roll_pitch_yaw(bot_matrix_to_quat(x0));
        rpy_guess = [r;p;y];
        tic
        [x,resnorm,residual,exitflag,output] = lsqnonlin(obj3,rpy_guess);
        used_time = toc;
        R_opt = bot_quat_to_matrix(bot_roll_pitch_yaw_to_quat(x));
        x1 = R_opt;
        cost = norm(obj3(x),2);

        if mode == "SB"
            A = zeros(3*(numTimes-1), numTimes+4);
        elseif mode == "S"
            A = zeros(3*(numTimes-1), numTimes+1);
        elseif mode == "B"
            A = zeros(3*(numTimes-1), numTimes+3);
        elseif mode == "Only"
            A = zeros(3*(numTimes-1), numTimes);
        end
    
        b = zeros(3*(numTimes-1), 1);
        for t = 1:numTimes-1
            if mode == "SB"
                A_t = zeros(3, numTimes+4);
            elseif mode == "S"
                A_t = zeros(3, numTimes+1);
            elseif mode == "B"
                A_t = zeros(3, numTimes+3);
            elseif mode == "Only"
                A_t = zeros(3, numTimes);
            end
            b_t = zeros(3,1);
        
            R_j0_jn0 = local_trajs(1).T(t).R;
            R_j0_jn1 = local_trajs(1).T(t+1).R;
            t_j0_jn0 = local_trajs(1).T(t).t;
            t_j0_jn1 = local_trajs(1).T(t+1).t;
        
            R_i0_in0 = local_trajs(2).T(t).R;
            R_i0_in1 = local_trajs(2).T(t+1).R;
            t_i0_in0 = local_trajs(2).T(t).t;
            t_i0_in1 = local_trajs(2).T(t+1).t;
        
            bi0 = bearings_data(1).bearing(t).x;
            bi1 = bearings_data(1).bearing(t+1).x;
        
            if mode == "SB"
                A_t(1:3, 1:3) = x1 * (R_i0_in0 - R_i0_in1);
                A_t(1:3, 4) = -x1*(t_i0_in1 - t_i0_in0);
                A_t(1:3, 4+t) = - R_j0_jn0 * bi0;
                A_t(1:3, 4+t+1) =  R_j0_jn1 * bi1;
                b_t = -(t_j0_jn1 - t_j0_jn0); 
            elseif mode == "S"
                A_t(1:3, 1) = -x1*(t_i0_in1 - t_i0_in0);
                A_t(1:3, 1+t) = - R_j0_jn0 * bi0;
                A_t(1:3, 1+t+1) =  R_j0_jn1 * bi1;
                b_t = -(t_j0_jn1 - t_j0_jn0); 
            elseif mode == "B"
                A_t(1:3, 1:3) = x1 * (R_i0_in0 - R_i0_in1);
                A_t(1:3, 3+t) = - R_j0_jn0 * bi0;
                A_t(1:3, 3+t+1) =  R_j0_jn1 * bi1;
                b_t = -(t_j0_jn1 - t_j0_jn0)+x1*(t_i0_in1 - t_i0_in0);
            elseif mode == "Only"
                A_t(1:3, t) = - R_j0_jn0 * bi0;
                A_t(1:3, t+1) =  R_j0_jn1 * bi1;
                b_t = -(t_j0_jn1 - t_j0_jn0)+x1*(t_i0_in1 - t_i0_in0);
            end
        
            A((t-1)*3+1 : 3*t, :) = A_t;
            b((t-1)*3+1 : 3*t, :) = b_t;
        end
    
        x2 = inv(A'*A)*A'*b;

        t_opt = zeros(3,1);
        for t = 1:numTimes
            R_j0_jn0 = local_trajs(1).T(t).R;
            t_j0_jn0 = local_trajs(1).T(t).t;
        
            R_i0_in0 = local_trajs(2).T(t).R;
            t_i0_in0 = local_trajs(2).T(t).t;
            
            bi0 = bearings_data(1).bearing(t).x;
            
            if mode == "SB"
                OP_bias = x2(1:3,:);
                OP_s = x2(4,:);
                d = x2(5:end,:);
            elseif mode == "S"
                OP_bias = zeros(3,1);
                OP_s = x2(1,:);
                d = x2(2:end,:);
            elseif mode == "B"
                OP_bias = x2(1:3,:);
                OP_s = 1.0;
                d = x2(4:end,:);
            elseif mode == "Only"
                OP_bias = zeros(3,1);
                OP_s = 1.0;
                d = x2;
            end

            depth = d(t);
            t_opt = t_opt + R_j0_jn0 * depth * bi0 + t_j0_jn0 - OP_s*x1 * (R_i0_in0 * OP_bias + t_i0_in0);
        end
        t_opt = t_opt / numTimes; 

        output.cost_table(i,j) = cost;
        output.times(i,j) = used_time;
        output.iteration(i,j) = output.iterations;
        output.l2error_R(i,j) = norm(R_opt-T_gt.R,2);
        output.l2error_t(i,j) = norm(t_opt-T_gt.t,2);
    end
end

end

function [cost] = costMarg(rpy, local_trajs, bearings_data, mode)
    R = bot_quat_to_matrix(bot_roll_pitch_yaw_to_quat(rpy));
    x1 = R;
    numTimes = size(local_trajs(1).poses,2);
    %% Given R* we can get Other*
    if mode == "SB"
        A = zeros(3*(numTimes-1), numTimes+4);
    elseif mode == "S"
        A = zeros(3*(numTimes-1), numTimes+1);
    elseif mode == "B"
        A = zeros(3*(numTimes-1), numTimes+3);
    elseif mode == "Only"
        A = zeros(3*(numTimes-1), numTimes);
    end

    b = zeros(3*(numTimes-1), 1);
    for t = 1:numTimes-1
        if mode == "SB"
            A_t = zeros(3, numTimes+4);
        elseif mode == "S"
            A_t = zeros(3, numTimes+1);
        elseif mode == "B"
            A_t = zeros(3, numTimes+3);
        elseif mode == "Only"
            A_t = zeros(3, numTimes);
        end
        b_t = zeros(3,1);
    
        R_j0_jn0 = local_trajs(1).T(t).R;
        R_j0_jn1 = local_trajs(1).T(t+1).R;
        t_j0_jn0 = local_trajs(1).T(t).t;
        t_j0_jn1 = local_trajs(1).T(t+1).t;
    
        R_i0_in0 = local_trajs(2).T(t).R;
        R_i0_in1 = local_trajs(2).T(t+1).R;
        t_i0_in0 = local_trajs(2).T(t).t;
        t_i0_in1 = local_trajs(2).T(t+1).t;
    
    
        bi0 = bearings_data(1).bearing(t).x;
        bi1 = bearings_data(1).bearing(t+1).x;
    
        if mode == "SB"
            A_t(1:3, 1:3) = x1 * (R_i0_in0 - R_i0_in1);
            A_t(1:3, 4) = -x1*(t_i0_in1 - t_i0_in0);
            A_t(1:3, 4+t) = - R_j0_jn0 * bi0;
            A_t(1:3, 4+t+1) =  R_j0_jn1 * bi1;
            b_t = -(t_j0_jn1 - t_j0_jn0); 
        elseif mode == "S"
            A_t(1:3, 1) = -x1*(t_i0_in1 - t_i0_in0);
            A_t(1:3, 1+t) = - R_j0_jn0 * bi0;
            A_t(1:3, 1+t+1) =  R_j0_jn1 * bi1;
            b_t = -(t_j0_jn1 - t_j0_jn0); 
        elseif mode == "B"
            A_t(1:3, 1:3) = x1 * (R_i0_in0 - R_i0_in1);
            A_t(1:3, 3+t) = - R_j0_jn0 * bi0;
            A_t(1:3, 3+t+1) =  R_j0_jn1 * bi1;
            b_t = -(t_j0_jn1 - t_j0_jn0)+x1*(t_i0_in1 - t_i0_in0);
        elseif mode == "Only"
            A_t(1:3, t) = - R_j0_jn0 * bi0;
            A_t(1:3, t+1) =  R_j0_jn1 * bi1;
            b_t = -(t_j0_jn1 - t_j0_jn0)+x1*(t_i0_in1 - t_i0_in0);
        end
    
        A((t-1)*3+1 : 3*t, :) = A_t;
        b((t-1)*3+1 : 3*t, :) = b_t;
        
    end

    x2 = inv(A'*A)*A'*b;
    
    cost = zeros(3,1);
    for t = 1:numTimes-1
    
        if mode == "SB"
            A_t = zeros(3, numTimes+4);
        elseif mode == "S"
            A_t = zeros(3, numTimes+1);
        elseif mode == "B"
            A_t = zeros(3, numTimes+3);
        elseif mode == "Only"
            A_t = zeros(3, numTimes);
        end
        b_t = zeros(3,1);
    
        R_j0_jn0 = local_trajs(1).T(t).R;
        R_j0_jn1 = local_trajs(1).T(t+1).R;
        t_j0_jn0 = local_trajs(1).T(t).t;
        t_j0_jn1 = local_trajs(1).T(t+1).t;
    
        R_i0_in0 = local_trajs(2).T(t).R;
        R_i0_in1 = local_trajs(2).T(t+1).R;
        t_i0_in0 = local_trajs(2).T(t).t;
        t_i0_in1 = local_trajs(2).T(t+1).t;
    
    
        bi0 = bearings_data(1).bearing(t).x;
        bi1 = bearings_data(1).bearing(t+1).x;
    
        if mode == "SB"
            A_t(1:3, 1:3) = x1 * (R_i0_in0 - R_i0_in1);
            A_t(1:3, 4) = -x1*(t_i0_in1 - t_i0_in0);
            A_t(1:3, 4+t) = - R_j0_jn0 * bi0;
            A_t(1:3, 4+t+1) =  R_j0_jn1 * bi1;
            b_t = -(t_j0_jn1 - t_j0_jn0); 
        elseif mode == "S"
            A_t(1:3, 1) = -x1*(t_i0_in1 - t_i0_in0);
            A_t(1:3, 1+t) = - R_j0_jn0 * bi0;
            A_t(1:3, 1+t+1) =  R_j0_jn1 * bi1;
            b_t = -(t_j0_jn1 - t_j0_jn0); 
        elseif mode == "B"
            A_t(1:3, 1:3) = x1 * (R_i0_in0 - R_i0_in1);
            A_t(1:3, 3+t) = - R_j0_jn0 * bi0;
            A_t(1:3, 3+t+1) =  R_j0_jn1 * bi1;
            b_t = -(t_j0_jn1 - t_j0_jn0)+x1*(t_i0_in1 - t_i0_in0);
        elseif mode == "Only"
            A_t(1:3, t) = - R_j0_jn0 * bi0;
            A_t(1:3, t+1) =  R_j0_jn1 * bi1;
            b_t = -(t_j0_jn1 - t_j0_jn0)+x1*(t_i0_in1 - t_i0_in0);
        end
    
        cost = cost + A_t * x2 - b_t;
        
    end

end