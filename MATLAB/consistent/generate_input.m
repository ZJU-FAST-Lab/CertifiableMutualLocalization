function [auxiliary_data,data_gt] = generate_input(local_trajs, bearings_data, world_trajs, param)

    numRobots = size(local_trajs,2);
    numTimes = size(local_trajs(1).poses,2);
    
    dim_BA = (3*numRobots)^2;
    dim_distance_variable = (numRobots-1)*numRobots*numTimes;
    all_variable_size = dim_BA + dim_distance_variable + 1;
    
    %% distance groundtruth
    distance_gt = zeros(dim_distance_variable, 1);
    for i = 1:numRobots
        for j = 1:numRobots
            if(i==j)
                continue
            end
            for t = 1:numTimes
                b_i_t1 = bearings_data(i,j).bearing(t).x;
                d_i_t1 = bearings_data(i,j).depth(t);
                if(j > i)
                    virtual_j = j-1;
                else
                    virtual_j = j;
                end
                index_distance_i_t1 = (i-1)*(numRobots-1)*numTimes + (virtual_j-1)*numTimes + t;
                distance_gt(index_distance_i_t1,1) = d_i_t1;
            end
        end
    end
    %% R,t groundtruth
    R_gt = zeros(3,3*numRobots);
    t_gt = zeros(3,numRobots);
    for i = 1:numRobots
        R_gt(:,3*(i-1)+1:3*(i-1)+3) = world_trajs(1).T(1).R' * world_trajs(i).T(1).R;
        t_gt(:, i) = world_trajs(1).T(1).R' * world_trajs(i).T(1).t; 
    end
    
    data_gt.distance_gt = distance_gt;
    data_gt.Z_gt = R_gt'*R_gt;
    data_gt.vecRTR_gt = vec(R_gt'*R_gt);
    data_gt.vect_gt = vec(t_gt);
    data_gt.all_gt = [data_gt.vecRTR_gt; 1; data_gt.distance_gt];
    data_gt.part_gt = [data_gt.vecRTR_gt; 1];
    data_gt.R_gt = R_gt(:,1:3)' * R_gt;
    data_gt.t_gt = t_gt;
    
    %% =============Construct Some Important Matrixes============
    M = zeros(all_variable_size, all_variable_size);
    BA_all = zeros(dim_BA, dim_BA);
    C_all = zeros(3*numRobots, 3*numRobots);
    D_all = 0;
    
    distance_b_static = zeros(dim_distance_variable,1);
    distance_b_wait = cell(dim_distance_variable,1);
    for i = 1:dim_distance_variable
        distance_b_wait{i,1} = zeros(3*numRobots, 3*numRobots);
    end
    distance_A = zeros(dim_distance_variable, dim_distance_variable);
    distance_A = sparse(distance_A);

    I = eye(3*numRobots);
    eps = 0;
    skip = 1;
    for i = 1:numRobots
        for j = 1:numRobots
            if(i==j)
                continue
            end
            for t = 1:numTimes - skip
                Q0 = zeros(3*numRobots, 3*numRobots);
                Ri  = world_trajs(i).T(1).R;  ti = world_trajs(i).T(1).t;
                Rj  = world_trajs(j).T(1).R;  tj = world_trajs(j).T(1).t;

                R_i_t1 = local_trajs(i).T(t).R;      t_i_t1 = local_trajs(i).T(t).t;
                R_i_t2 = local_trajs(i).T(t+skip).R; t_i_t2 = local_trajs(i).T(t+skip).t;

                R_j_t1 = local_trajs(j).T(t).R;      t_j_t1 = local_trajs(j).T(t).t;
                R_j_t2 = local_trajs(j).T(t+skip).R; t_j_t2 = local_trajs(j).T(t+skip).t;
    
                b_i_t1 = bearings_data(i,j).bearing(t).x;      d_i_t1 = bearings_data(i,j).depth(t);
                b_i_t2 = bearings_data(i,j).bearing(t+skip).x; d_i_t2 = bearings_data(i,j).depth(t+skip);
                
                t_i_t1_t2 = t_i_t2 - t_i_t1; t_j_t1_t2 = t_j_t2 - t_j_t1;
                g_i_t1 = R_i_t1 * b_i_t1; g_i_t2 = R_i_t2 * b_i_t2;
                
                % Choose Matrix
                C1 = zeros(3*numRobots,3);
                C2 = zeros(3*numRobots,3);
                C1(3*(i-1)+1:3*(i-1)+3,:) = eye(3);
                C2(3*(j-1)+1:3*(j-1)+3,:) = eye(3);

                
                if(param.mode == "Marginalization")
                    % Version: Marg
                    A_var =C1'*kron((C2*t_j_t1_t2)', I);
                    B_var = -t_i_t1_t2;
                    C_var = zeros(3,d_variable_size);
                    if(j > i)
                        virtual_j = j-1;
                    else
                        virtual_j = j;
                    end
                    index_distance_i_t1 = (i-1)*(numRobots-1)*numTimes + (virtual_j-1)*numTimes + t;
                    index_distance_i_t2 = index_distance_i_t1 + 1;
                    C_var(:,index_distance_i_t1) = R_i_t1 * b_i_t1;
                    C_var(:,index_distance_i_t2) = -R_i_t2 * b_i_t2;
                    all_var = [A_var B_var C_var];
                    M = M + all_var'*all_var;
                else
                    % Version: No Marg + Information
                    k_i_t1_t2 = cross(g_i_t1,g_i_t2);
                    K_12 = k_i_t1_t2 * k_i_t1_t2';

                    A_part = C2*(t_j_t1_t2*t_j_t1_t2')*C2' + eps * eye(3*numRobots);
                    B_part = C1*(k_i_t1_t2 * k_i_t1_t2')*C1' + eps * eye(3*numRobots);
                    C_part = -2*C2*(t_j_t1_t2*t_i_t1_t2')*K_12*C1';
                    D_part = t_i_t1_t2'*K_12*t_i_t1_t2;

                    BA_part = kron(B_part, A_part);
                    BA_all = BA_all + BA_part;
                    C_all = C_all + C_part;
                    D_all = D_all + D_part;
                end

                % for distance recovery
                if(j > i)
                    virtual_j = j-1;
                else
                    virtual_j = j;
                end
                index_distance_i_t1 = (i-1)*(numRobots-1)*numTimes + (virtual_j-1)*numTimes + t;
                index_distance_i_t2 = index_distance_i_t1 + skip;
                
                distance_A(index_distance_i_t1, index_distance_i_t1) = distance_A(index_distance_i_t1, index_distance_i_t1) + norm(g_i_t1);
                distance_A(index_distance_i_t1, index_distance_i_t2) = distance_A(index_distance_i_t1, index_distance_i_t2) + -g_i_t1' * g_i_t2;
                distance_A(index_distance_i_t2, index_distance_i_t1) = distance_A(index_distance_i_t2, index_distance_i_t1) + -g_i_t1' * g_i_t2;
                distance_A(index_distance_i_t2, index_distance_i_t2) = distance_A(index_distance_i_t2, index_distance_i_t2) + norm(g_i_t2);

                distance_b_wait{index_distance_i_t1,1} =  distance_b_wait{index_distance_i_t1,1} + -C2*t_j_t1_t2*(-g_i_t1)'*C1';
                distance_b_wait{index_distance_i_t2,1} =  distance_b_wait{index_distance_i_t2,1} + -C2*t_j_t1_t2*( g_i_t2)'*C1';

                distance_b_static(index_distance_i_t1,1) =  distance_b_static(index_distance_i_t1,1) + (-g_i_t1'*t_i_t1_t2);
                distance_b_static(index_distance_i_t2,1) =  distance_b_static(index_distance_i_t2,1) + g_i_t2'*t_i_t1_t2;
            end
        end
    end

    %% Marginalization
     if(param.mode == "Marginalization")
        M_rr11 = M(1:BA_dim+1, 1:BA_dim+1);
        M_dd = M(BA_dim+2:all_variable_size, BA_dim+2:all_variable_size);
        M_r1_d = M(1:BA_dim+1, BA_dim+2:all_variable_size);

        M_marg = M_rr11 - M_r1_d * pinv(M_dd) * M_r1_d';

        BA_all = M_marg(1:BA_dim,1:BA_dim);
        C_all = 2*vec_inv(M_marg(1:BA_dim,BA_dim+1));
        D_all = M_marg(end,end);
     end

    %% Set auxiliary data
    auxiliary_data.BA = BA_all;
    auxiliary_data.C = C_all;
    auxiliary_data.D = D_all;
    auxiliary_data.d = 3;
    auxiliary_data.n = numRobots;
    auxiliary_data.dim_distance_variable = dim_distance_variable;
    auxiliary_data.distance_A = distance_A;
    auxiliary_data.distance_b_wait = distance_b_wait;
    auxiliary_data.distance_b_static = distance_b_static;
end

