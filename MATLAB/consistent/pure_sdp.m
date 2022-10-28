function [result] = pure_sdp(data)
 %% Iteration
    numRobots = data.numRobots;
    auxiliary_data = data.auxiliary_data;
    dims = 3*numRobots;
    mark_tic = tic();

    
    % Set CVX options
    cvx_precision('high') 
    % cvx_solver sedumi
    cvx_solver SDPT3
    cvx_quiet(false) % necessary false to dump results
    dumpfile = sprintf('temp_cvxDump_%s',num2str(randi(1e10)));
    cvx_solver_settings( 'dumpfile', dumpfile ) % dump statistics for parsing
    cvx_begin SDP % Use CVX's SDP mode
    variable Z(dims, dims) symmetric
    variable t
    Z == semidefinite(dims);
    for i= 1:numRobots
    %     Z(i,i) == 1;
        Z(3*(i-1)+1,3*(i-1)+1) == 1;
        Z(3*(i-1)+2,3*(i-1)+2) == 1;
        Z(3*(i-1)+3,3*(i-1)+3) == 1;
    
        Z(3*(i-1)+1,3*(i-1)+2) == 0;
        Z(3*(i-1)+1,3*(i-1)+3) == 0;
    
        Z(3*(i-1)+2,3*(i-1)+1) == 0;
        Z(3*(i-1)+2,3*(i-1)+3) == 0;
    
        Z(3*(i-1)+3,3*(i-1)+1) == 0;
        Z(3*(i-1)+3,3*(i-1)+2) == 0;
    end
    cost = evaluate_objective_WYJ_by_X(Z, auxiliary_data);
    cost <= t
    tic()
    minimize(t)
    cvx_end
    t_sdp = toc();
    delete([dumpfile,'.mat'])
    cvx_solver_settings -clear dumpfile
    mark_toc = toc(mark_tic);
    [U,S,V] = svd(Z);
    sqrt_S = sqrt(S);
    Yopt  = sqrt_S(1:3,1:3) * V(:,1:3)';
    if(det(Yopt(:,1:3)) < 0)
        Yopt = -Yopt;
    end
    R_opt_SDP = zeros(3, 3*numRobots);
    for k = 1:numRobots
        Y = Yopt(:,3*(k-1)+1:3*(k-1)+3);
        [U,S,V] = svd(Y);
        R_opt_SDP(:, 3*(k-1) + 1 : 3*(k-1) + 3) = U*V'*det(U*V');
    end
    f_X = evaluate_objective_WYJ_by_X(Z, auxiliary_data);
    g_RTR = evaluate_objective_WYJ(R_opt_SDP, auxiliary_data);

    result.R_opt = R_opt_SDP(:,1:3)'*R_opt_SDP;
    result.f_X = f_X;
    result.g_RTR = g_RTR;
    result.refuse_error = norm(Yopt'*Yopt - Z);
    result.rounding_error =  norm(R_opt_SDP - Yopt);
    result.time = mark_toc;
    result.l2error = norm(result.R_opt - data.data_gt.R_gt);
    result.X = Z;
end