function [output] = benchmark_tro(local_trajs, bearings_data, data_gt, tro_param, gt_T, P , mode, resolution, if_opt_with_gt)


T_gt = data_gt.T(1);
R_gt = T_gt.R;
resolution = tro_param.resolution;
mode = tro_param.mode;
opt_using_gt = tro_param.opt_using_gt;

output.cost_table = zeros(resolution, resolution);
output.times = zeros(resolution, resolution);
output.iteration = zeros(resolution, resolution);
output.l2error_R = zeros(resolution, resolution);
output.l2error_t = zeros(resolution, resolution);

for i = 1:resolution
    roll = i * 2 * pi / resolution - pi;
    for j = 1:resolution
        [i,j]
        pitch = j * 2 * pi/ resolution- pi;
        rpy = [roll;pitch;0];
        disturb = bot_quat_to_matrix(bot_roll_pitch_yaw_to_quat(rpy));
        x1  = disturb * R_gt; 
        
        if(resolution == 1)
            x1 = rnd.rot();
        end
        if(opt_using_gt)
            x1 = R_gt;
        end

        numTimes = size(local_trajs(1).poses,2);
        cost = 1e10;
        tic
        for iter = 1:10000
           %% fix x1 
            if mode == "SB"
                A = zeros(3*(numTimes-1), numTimes+4);
            elseif  mode == "S"
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
        
           %% fix x2
            if mode == "SB"
                P = x2(1:3,:);
                s = x2(4,:);
            elseif mode == "S"
                s = x2(1,:);
            elseif mode == "B"
                P = x2(1:3,:);
            elseif mode == "Only"
            end
        
            C = zeros(3, numTimes-1);
            D = zeros(3, numTimes-1);
            for t = 1:numTimes-1
                if mode == "SB"
                    Z_t_0 = x2(4+t,:);
                    Z_t_1 = x2(4+t+1,:);
                elseif mode == "S"
                    Z_t_0 = x2(1+t,:);
                    Z_t_1 = x2(1+t+1,:);
                elseif mode == "B"
                    Z_t_0 = x2(3+t,:);
                    Z_t_1 = x2(3+t+1,:);
                elseif mode == "Only"
                    Z_t_0 = x2(t,:);
                    Z_t_1 = x2(t+1,:);
                end
        
                C_t = zeros(3,1);
                D_t = zeros(3,1);
                
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
                    C_t = (R_i0_in0 - R_i0_in1) * P - s * (t_i0_in1 - t_i0_in0);
                    D_t = R_j0_jn0 * Z_t_0 * bi0 - R_j0_jn1 * Z_t_1 * bi1 - (t_j0_jn1 - t_j0_jn0);
                elseif mode == "S"
                    C_t = - s * (t_i0_in1 - t_i0_in0);
                    D_t = R_j0_jn0 * Z_t_0 * bi0 - R_j0_jn1 * Z_t_1 * bi1 - (t_j0_jn1 - t_j0_jn0);
                elseif mode == "B"
                    C_t = (R_i0_in0 - R_i0_in1) * P -  (t_i0_in1 - t_i0_in0);
                    D_t = R_j0_jn0 * Z_t_0 * bi0 - R_j0_jn1 * Z_t_1 * bi1 - (t_j0_jn1 - t_j0_jn0);
                elseif mode == "Only"
                    C_t = - (t_i0_in1 - t_i0_in0);
                    D_t = R_j0_jn0 * Z_t_0 * bi0 - R_j0_jn1 * Z_t_1 * bi1 - (t_j0_jn1 - t_j0_jn0);
                end
        
                C(:,t) = C_t; 
                D(:,t) = D_t;
            end
            
            [U,S,V] = svd(C*D');
            G = eye(3);
            G(3,3) = det(U*V);
            x1 = V*G*U';
        
            cost1 = norm(x1 * C - D,2);
            if(abs(cost1 - cost) / cost1 < 0.01)
                cost = cost1;
                break
            end
            cost = cost1;
        end
        
        t_opt = zeros(3,1);
        for t = 1:numTimes
            R_j0_jn0 = local_trajs(1).T(t).R;
            t_j0_jn0 = local_trajs(1).T(t).t;
            bi0 = bearings_data(1).bearing(t).x;
        
            R_i0_in0 = local_trajs(2).T(t).R;
            t_i0_in0 = local_trajs(2).T(t).t;
            
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
        R_opt = x1;
        
        output.cost_table(i,j) = cost;
        output.times(i,j) = toc;
        output.iteration(i,j) = iter;
        output.l2error_R(i,j) = norm(R_opt-T_gt.R,2);
        output.l2error_t(i,j) = norm(t_opt-T_gt.t,2);
    end
end
end