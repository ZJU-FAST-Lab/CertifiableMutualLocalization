function [output] = simple_method(local_trajs, bearings_data, data_gt, param)

numRobots = size(local_trajs,2);
numTimes = size(local_trajs(1).poses,2);

order_gt = data_gt.order;
T_gt = data_gt.T;
s_gt = data_gt.s;

%% Some definition
% s: scale
% p1, p2, p3:three elements of bias 
% k: theta_ij, i for traj, j for bearings
weight = 1;
% dimension of variables, the variables includes 
% k*s*R, k*p1*R, k*p2*R, k*p3*R [k,s1,s2,s3 for 4, R for 9, k for (numRobots-1)]
% and all distance variables (for numTimes)
dim_for_one_robot = 9*4*(numRobots-1) + numTimes; 
dim = (numRobots-1)*dim_for_one_robot + 1; %% +1 is for homogeneous form

dim_nD_for_one_robot = 9*4*(numRobots-1); % eliminate distance variables
dim_nD = (numRobots-1)*dim_nD_for_one_robot + 1;

dim_lift_for_one_robot = 9*4*(numRobots-1 + 1); % add additional variables: s*R p1*R p2*R p3*R
dim_lift = (numRobots-1)*dim_lift_for_one_robot + 1;

dim_lift_plus_K = dim_lift + (numRobots-1)^2; % add data correspondence variables: Theta
dim_lift_plus_K_SP = dim_lift + (numRobots-1)^2 + (numRobots-1)*4; %add variables:s, p1, p2, p3
dim_lift_plus_K_SP_KSP = dim_lift + (numRobots-1)^2 + (numRobots-1)*4 + (numRobots-1)^2*4;%add variables: k*s k*p1 k*p2 k*p3

Q = zeros(dim,dim); % for QCQP: min x'Qx
gd_x = zeros(dim,1); % groudtruth state
gd_x_lift = zeros(dim_lift_plus_K_SP_KSP,1); % lifted groudtruth state

%permutation matrix
K = zeros(numRobots-1, numRobots-1);
for i = 1:numRobots-1
    K(order_gt(1,i),i) = 1;
end
K = K';

%% Construct groudtruth state gd_x and gd_x_lift
%the index of the variable 1
index_one = dim_nD;
gd_x(index_one) = 1;
gd_x_lift(index_one) = 1;
for r = 2:numRobots
    for k = 2:numRobots
        % indexes of variables
        index_kvecsR = (r-2)*4*9*(numRobots-1)+(k-2)*4*9;
        index_kvecp1R = (r-2)*4*9*(numRobots-1)+(k-2)*4*9+9;
        index_kvecp2R = (r-2)*4*9*(numRobots-1)+(k-2)*4*9+18;
        index_kvecp3R = (r-2)*4*9*(numRobots-1)+(k-2)*4*9+27;
        
        index_vecsR = dim_nD+(k-2)*4*9;
        index_vecp1R = dim_nD+(k-2)*4*9+9;
        index_vecp2R = dim_nD+(k-2)*4*9+18;
        index_vecp3R = dim_nD+(k-2)*4*9+27;

        index_k = dim_lift+(r-2)*(numRobots-1)+(k-1); %k(r)(k)
        
        index_s = dim_lift_plus_K + (k-2)*4+1; %s(k)
        index_p1 = dim_lift_plus_K + (k-2)*4+2; %p1(k)
        index_p2 = dim_lift_plus_K + (k-2)*4+3; %p2(k)
        index_p3 = dim_lift_plus_K + (k-2)*4+4; %p3(k)

        index_ks = dim_lift_plus_K_SP + (r-2)*(numRobots-1)*4+(k-2)*4+1; %k(r)(k) * s(k)
        index_kp1 = dim_lift_plus_K_SP + (r-2)*(numRobots-1)*4+(k-2)*4+2; %k(r)(k) * p1(k)
        index_kp2 = dim_lift_plus_K_SP + (r-2)*(numRobots-1)*4+(k-2)*4+3; %k(r)(k) * p2(k)
        index_kp3 = dim_lift_plus_K_SP + (r-2)*(numRobots-1)*4+(k-2)*4+4; %k(r)(k) * p3(k)
        
       % groundtruth value
        ks = s_gt(k-1);
        R_gt = T_gt(k-1).R;
        gt_bias = local_trajs(k).bias;
        p1 = gt_bias(1); p2 = gt_bias(2); p3 = gt_bias(3);
        
       % constrict gd_x
        gd_x(index_kvecsR+1:index_kvecsR+9,:) = vec(ks*R_gt)*K(r-1,k-1);
        gd_x(index_kvecp1R+1:index_kvecp1R+9,:) = vec(p1*R_gt)*K(r-1,k-1);
        gd_x(index_kvecp2R+1:index_kvecp2R+9,:) = vec(p2*R_gt)*K(r-1,k-1);
        gd_x(index_kvecp3R+1:index_kvecp3R+9,:) = vec(p3*R_gt)*K(r-1,k-1);

       % construct gd_x_lift
        gd_x_lift(index_kvecsR+1:index_kvecsR+9,:) = vec(ks*R_gt)*K(r-1,k-1);
        gd_x_lift(index_kvecp1R+1:index_kvecp1R+9,:) = vec(p1*R_gt)*K(r-1,k-1);
        gd_x_lift(index_kvecp2R+1:index_kvecp2R+9,:) = vec(p2*R_gt)*K(r-1,k-1);
        gd_x_lift(index_kvecp3R+1:index_kvecp3R+9,:) = vec(p3*R_gt)*K(r-1,k-1);
        
        gd_x_lift(index_vecsR+1:index_vecsR+9,:) = vec(ks*R_gt);
        gd_x_lift(index_vecp1R+1:index_vecp1R+9,:) = vec(p1*R_gt);
        gd_x_lift(index_vecp2R+1:index_vecp2R+9,:) = vec(p2*R_gt);
        gd_x_lift(index_vecp3R+1:index_vecp3R+9,:) = vec(p3*R_gt);
        
        gd_x_lift(index_k,:) = K(r-1,k-1);
        
        gd_x_lift(index_s,:) = ks;
        gd_x_lift(index_p1,:) = p1;
        gd_x_lift(index_p2,:) = p2;
        gd_x_lift(index_p3,:) = p3;
        
        gd_x_lift(index_ks,:) = ks*K(r-1,k-1);
        gd_x_lift(index_kp1,:) = p1*K(r-1,k-1);
        gd_x_lift(index_kp2,:) = p2*K(r-1,k-1);
        gd_x_lift(index_kp3,:) = p3*K(r-1,k-1);
    end
   
    % constrict gd_x
    for i = 1:numTimes
        index_distance = dim_nD + numTimes*(r-2)+i;
        d = bearings_data(r-1).depth(i);
        gd_x(index_distance,:) = d;
    end
end


%% Construct the cost matrix Q for QCQP: min x'Qx
for r = 2:numRobots
    for skip = 1:1 % ||frame_i - frame_j||
        for i = 1:numTimes-skip
            vec1 = zeros(3,dim);
            R_j0_jn0 = local_trajs(1).T(i).R; t_j0_jn0 = local_trajs(1).T(i).t;
            R_j0_jn1 = local_trajs(1).T(i+skip).R; t_j0_jn1 = local_trajs(1).T(i+skip).t;
            bi0 = bearings_data(r-1).bearing(i).x; bi1 = bearings_data(r-1).bearing(i+skip).x;
            Rn0b0 = R_j0_jn0*bi0;  Rn1b1 = R_j0_jn1*bi1;
            tA = t_j0_jn1 - t_j0_jn0;
            
            for k = 2:numRobots
                R_i0_in0 = local_trajs(k).T(i).R; t_i0_in0 = local_trajs(k).T(i).t;
                R_i0_in1 = local_trajs(k).T(i+skip).R; t_i0_in1 = local_trajs(k).T(i+skip).t;
                v_lift_0 = [t_i0_in0;vec(R_i0_in0)]; v_lift_1 = [t_i0_in1;vec(R_i0_in1)];

                index_kvecsR = (r-2)*4*9*(numRobots-1)+(k-2)*4*9;
                vec1(:,index_kvecsR+1:index_kvecsR+36) = kron((v_lift_1-v_lift_0)' ,eye(3));
            end
            
            index_distance = dim_nD + numTimes*(r-2)+i;
            index_distance2 = dim_nD + numTimes*(r-2)+i+skip;
            vec1(:,index_one) = -tA;
            vec1(:,index_distance) = Rn0b0;
            vec1(:,index_distance2) = -Rn1b1;

            N = vec1'*weight*vec1;
            Q = N + Q;
        end
    end
end

%% Marginalization
Q_nZ_nZ = Q(1:dim_nD, 1:dim_nD);
Q_nZ_Z = Q(1:dim_nD, dim_nD+1:end);
Q_Z_nZ = Q(dim_nD+1:end, 1:dim_nD);
Q_Z_Z = Q(dim_nD+1:end, dim_nD+1:end);
Q_marg = Q_nZ_nZ - Q_nZ_Z * pinv(Q_Z_Z) * Q_Z_nZ;

%% Add additional variable and constrains
Q_lift = zeros(dim_lift_plus_K_SP_KSP,dim_lift_plus_K_SP_KSP);
Q_lift(1:dim_nD, 1:dim_nD) = Q_marg;

%% =======Begin SDP=======
% Set CVX options
cvx_precision('best') 
% cvx_solver sedumi
% cvx_solver SDPT3
cvx_quiet(false) % necessary false to dump results
dumpfile = sprintf('temp_cvxDump_%s',num2str(randi(1e10)));
cvx_solver_settings( 'dumpfile', dumpfile ) % dump statistics for parsing

cvx_begin SDP % Use CVX's SDP mode
variable Z(dim_lift_plus_K_SP_KSP, dim_lift_plus_K_SP_KSP) symmetric

Z == semidefinite(dim_lift_plus_K_SP_KSP);
dual variable lambda{3000}
dual_num = 1;

%% =======Add constrains=======
%% k*k = k*1
% k:theta_ij
for r = 2:numRobots
    for k = 2:numRobots
        index_k = dim_lift+(r-2)*(numRobots-1)+(k-1); %k(r)(k)

        P = zeros(dim_lift_plus_K_SP_KSP,dim_lift_plus_K_SP_KSP);
        P(index_k, index_k) = 1;
        P(index_k, index_one) = -1;
        P = symmetrize(P);
        
        trace(P*Z) == 0:lambda{dual_num}; 
        dual_num = dual_num + 1;
    end
end
disp("k*k = k*1")
dual_num

%% SUM(K)=1
%% k1+k2+...+kn=1 row
%% k1+k2+...+kn=1 col
% k:theta_ij
for r = 2:numRobots
    % ROW!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! r first and then k
    P = zeros(dim_lift_plus_K_SP_KSP,dim_lift_plus_K_SP_KSP);
    for k  = 2:numRobots
        index_k = dim_lift+(r-2)*(numRobots-1)+(k-1); %k(r)(k)
        P(index_one, index_k) = 1;
    end
    P(index_one, index_one) = -1;
    P = symmetrize(P);
    
    trace(P*Z) == 0:lambda{dual_num}; 
    dual_num = dual_num + 1;

    % COL!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! k first and then r
    P = zeros(dim_lift_plus_K_SP_KSP,dim_lift_plus_K_SP_KSP);
    for k  = 2:numRobots
        index_k = dim_lift+(k-2)*(numRobots-1)+(r-1);
        P(index_one, index_k) = 1;
    end
    P(index_one, index_one) = -1;
    P = symmetrize(P);

    trace(P*Z) == 0:lambda{dual_num}; 
    dual_num = dual_num + 1;
end
disp("k1+k2+...+kn=1")
dual_num


%% Rotation Constraints for vec(xR) 
% x=s,p1,p2,p3
for r = 2:numRobots
    for i = 1:4 %i = 1:s | 2:p1 | 3:p2 | 4:p3
        index_vecxR = dim_nD+(r-2)*4*9 + (i-1)*9;
        index_x = dim_lift_plus_K + (r-2)*4+i;
        % (xR)^T (xR) = x^2 I
        for m = 1:3
            for n=m:3
                P = zeros(dim_lift_plus_K_SP_KSP,dim_lift_plus_K_SP_KSP);
                P(index_vecxR+1:index_vecxR+9, index_vecxR+1:index_vecxR+9) = kron(symmetrize(canvec33(m,n)),eye(3));
                P(index_x, index_x) = -(m==n);
                P = symmetrize(P);
                
                trace(P*Z) == 0:lambda{dual_num}; 
                dual_num = dual_num + 1;
            end
        end

        % (xR) (xR)^T = x^2 I
        for m = 1:3
            for n=m:3
                P = zeros(dim_lift_plus_K_SP_KSP,dim_lift_plus_K_SP_KSP);
                P(index_vecxR+1:index_vecxR+9, index_vecxR+1:index_vecxR+9) = kron(eye(3), symmetrize(canvec33(m,n)));
                P(index_x, index_x) = -(m==n);
                P = symmetrize(P);
                
                trace(P*Z) == 0:lambda{dual_num}; 
                dual_num = dual_num + 1;
            end
        end
        
        % (xR)(:,i) x (xR)(:,j) = x * (xR)(:,k), (i,j,k)=(1,2,3),(2,3,1)or(3,2,1)
        for r_i = 1:3
            r_j = mod(r_i,3)+1;
            r_k = mod(r_i+1,3)+1;
            for r_a = 1:3
                P = zeros(dim_lift_plus_K_SP_KSP,dim_lift_plus_K_SP_KSP);
                P(index_vecxR+1:index_vecxR+9, index_vecxR+1:index_vecxR+9) = kron(skew_symmetrize(canvec33(r_i,r_j)),skew(canvec3(r_a)));
                P(index_x, index_vecxR+1:index_vecxR+9) = 0.5 * kron(canvec3(r_k),canvec3(r_a))';
                P(index_vecxR+1:index_vecxR+9, index_x) = 0.5 * kron(canvec3(r_k),canvec3(r_a));
                P = symmetrize(P);
                
                trace(P*Z) == 0:lambda{dual_num}; 
                dual_num = dual_num + 1;
            end
        end
    end
end
disp("Rotation Constraints for vecxR")
dual_num

%% Rotation for k*vec(xR) 
% k:theta_ij
% x:s,p1,p2,p3
for r = 2:numRobots
    for k = 2:numRobots
        for i = 1:4 %i = 1:s | 2:p1 | 3:p2 | 4:p3
            index_kvecxR = (r-2)*4*9*(numRobots-1)+(k-2)*4*9 + (i-1)*9;
            index_kx = dim_lift_plus_K_SP + (r-2)*(numRobots-1)*4+(k-2)*4+i; %k(r)(k) * s(k)
            % (kxR)^T (kxR) = (kx)^2 I
            for m = 1:3
                for n=m:3
                    P = zeros(dim_lift_plus_K_SP_KSP,dim_lift_plus_K_SP_KSP);
                    P(index_kvecxR+1:index_kvecxR+9, index_kvecxR+1:index_kvecxR+9) = kron(symmetrize(canvec33(m,n)),eye(3));
                    P(index_kx, index_kx) = -(m==n);
                    P = symmetrize(P);
                    
                    trace(P*Z) == 0:lambda{dual_num}; 
                    dual_num = dual_num + 1;
                end
            end
    
            % (kxR) (kxR)^T = (kx)^2 I
            for m = 1:3
                for n=m:3
                    P = zeros(dim_lift_plus_K_SP_KSP,dim_lift_plus_K_SP_KSP);
                    P(index_kvecxR+1:index_kvecxR+9, index_kvecxR+1:index_kvecxR+9) = kron(eye(3), symmetrize(canvec33(m,n)));
                    P(index_kx, index_kx) = -(m==n);
                    P = symmetrize(P);
                    
                    trace(P*Z) == 0:lambda{dual_num}; 
                    dual_num = dual_num + 1;
                end
            end
            
            % (kxR)(:,i) x (kxR)(:,j) = kx * (kxR)(:,k), (i,j,k)=(1,2,3),(2,3,1)or(3,2,1)
            for r_i = 1:3
                r_j = mod(r_i,3)+1;
                r_k = mod(r_i+1,3)+1;
                for r_a = 1:3
                    P = zeros(dim_lift_plus_K_SP_KSP,dim_lift_plus_K_SP_KSP);
                    P(index_kvecxR+1:index_kvecxR+9, index_kvecxR+1:index_kvecxR+9) = kron(skew_symmetrize(canvec33(r_i,r_j)),skew(canvec3(r_a)));
                    P(index_kx, index_kvecxR+1:index_kvecxR+9) = 0.5 * kron(canvec3(r_k),canvec3(r_a))';
                    P(index_kvecxR+1:index_kvecxR+9, index_kx) = 0.5 * kron(canvec3(r_k),canvec3(r_a));
                    P = symmetrize(P);
                    
                    trace(P*Z) == 0:lambda{dual_num}; 
                    dual_num = dual_num + 1;
                end
            end
        end
    end
end
disp("Rotation Constraints for kvecxR")
dual_num

%% Some Equality Constraints

%% kx = k*x
% k:theta_ij
% x:s,p1,p2,p3
for r = 2:numRobots
    for k = 2:numRobots
        index_k = dim_lift+(r-2)*(numRobots-1)+(k-1); %k(r)(k)
        for i = 1:4 %i = 1:s | 2:p1 | 3:p2 | 4:p3
            index_x = dim_lift_plus_K + (k-2)*4 + i;
            index_kx = dim_lift_plus_K_SP + (r-2)*(numRobots-1)*4 + (k-2)*4 + i; %k(r)(k) * x(k)

            P = zeros(dim_lift_plus_K_SP_KSP,dim_lift_plus_K_SP_KSP);
            P(index_kx, index_one) = 1;
            P(index_k, index_x) = -1;           
            P = symmetrize(P);

            trace(P*Z) == 0:lambda{dual_num}; 
            dual_num = dual_num + 1;
        end
    end
end
disp("k*x = kx")
dual_num

%% kvecxR = k*vecxR
% k:theta_ij
% x:s,p1,p2,p3
for r = 2:numRobots
    for k = 2:numRobots
        for i = 1:4 %i = 1:s | 2:p1 | 3:p2 | 4:p3
            index_kvecxR = (r-2)*4*9*(numRobots-1)+(k-2)*4*9 + (i-1)*9;
            index_vecxR = dim_nD+(k-2)*4*9 + (i-1)*9;
            index_k = dim_lift+(r-2)*(numRobots-1)+(k-1);

            for j = 1:9 %vec(xR)
                P = zeros(dim_lift_plus_K_SP_KSP,dim_lift_plus_K_SP_KSP);
                P(index_k, index_vecxR+j) = 1;
                P(index_one, index_kvecxR+j) = -1;
                P = symmetrize(P);
                
                trace(P*Z) == 0:lambda{dual_num}; 
                dual_num = dual_num + 1;
            end
        end
        
    end
end
disp("kvecxR = k*vecxR")
dual_num

%% y=1 constraint
Z(index_one,index_one) == 1:lambda{dual_num}; 

%% Optimize !
minimize(trace(Q_lift*Z))
tic
cvx_end
toc;

delete([dumpfile,'.mat'])
cvx_solver_settings -clear dumpfile

%% ======PostProcess======

[U,S,V] = svd(Z(1:dim_nD, 1:dim_nD));
S = diag(S);

z_star = sqrt(S(1)) * U(:,1);
z_star = z_star / z_star(index_one);

%1. get permutation matrix K_opt based on s
K_opt_float = zeros(numRobots-1, numRobots-1);
for r = 2:numRobots
    for k = 2:numRobots
         index_kvecsR = (r-2)*4*9*(numRobots-1)+(k-2)*4*9;
         kvecsR = z_star(index_kvecsR+1:index_kvecsR+9,:);
         ksR = reshape(kvecsR,3,3);
         ks = norm(ksR(:,1),2);
    
         K_opt_float(r-1,k-1) = ks;
    end
end

K_opt = zeros(numRobots-1, numRobots-1);
for r = 2:numRobots
    max = -1e10;
    max_index = -1;
    for k = 2:numRobots 
        if(K_opt(r-1,k-1) > max)
            max = K_opt_float(r-1,k-1);
            max_index = k;
        end
    end
    K_opt(r-1,max_index-1) = 1;
end

%2. get optimized depth
depth_opt = -inv(Q_Z_Z) * Q_Z_nZ * z_star;

%3. get optimized s,R,and bias and calculate error
bias_opt = zeros(3,numRobots-1);
s_opt = zeros(1,numRobots-1);
t_opt = zeros(3,numRobots-1);
R_opt = zeros(3,3*(numRobots-1));
l2error.R = 0; l2error.t = 0;
for i = 2:numRobots
    index = (find(K_opt(i-1,:)==1));
    index_kvecsR = (i-2)*4*9*(numRobots-1)+(index-1)*4*9;
    kvecsR = z_star(index_kvecsR+1:index_kvecsR+9,:);
    ksR = reshape(kvecsR,3,3);
    ks = norm(ksR(:,1),2);
    s_opt(1,index) = ks;
    R_opt_1 = ksR/ks;
    R_opt(:,3*(index-1)+1:3*(index-1)+3) = R_opt_1;
    
    for j = 1:3
         kvecpR = z_star(index_kvecsR+j*9+1:index_kvecsR+j*9+9,:);
         kpR = reshape(kvecpR,3,3);
         p = det(kpR);
         bias_opt(j,index) = norm(kpR(:,1),2) * p / abs(p);
    end

    t_opt_1 = zeros(3,1);
    for t = 1:numTimes
        R_j0_jn0 = local_trajs(1).T(t).R; t_j0_jn0 = local_trajs(1).T(t).t;
        R_i0_in0 = local_trajs(index+1).T(t).R; t_i0_in0 = local_trajs(index+1).T(t).t;
        bi0 = bearings_data(i-1).bearing(t).x;
        depth = depth_opt((index-1)*numTimes+t);
        t_opt_1 = t_opt_1 + R_j0_jn0 * depth * bi0 + t_j0_jn0 - ks * R_opt_1 * (R_i0_in0 * bias_opt(:,index) / ks + t_i0_in0);
    end
    t_opt_1 = t_opt_1 / numTimes; 
    t_opt(:,index) = t_opt_1;

    l2error.R = l2error.R + norm(R_opt_1 - T_gt(index).R,2);
    l2error.t = l2error.t + norm(t_opt_1 - T_gt(index).t,2);
end
l2error.R = l2error.R / (numRobots - 1);
l2error.t = l2error.t / (numRobots - 1);

%% Pack Result
output.error = l2error;
output.K_opt = K_opt;
output.R_opt = R_opt;
output.t_opt = t_opt;
output.bias_opt = bias_opt;
output.s_opt = s_opt;
output.depth_opt = depth_opt;
end