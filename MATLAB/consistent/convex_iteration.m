function [result] = convex_iteration(data)
 %% Iteration
    numRobots = data.numRobots;
    auxiliary_data = data.auxiliary_data;
    tolerance = data.tolerance;

    dims = 3*numRobots;
    C = eye(dims);
    scale = 100000;

    record = [];
    max_iter = 100;
    mark_tic = tic();
    for iter = 1:max_iter
        disp(sprintf('============%d==================)', iter)); 
        cvx_precision('high') 
        cvx_quiet(false) % necessary false to dump results
        dumpfile = sprintf('temp_cvxDump_%s',num2str(randi(1e10)));
        cvx_solver_settings( 'dumpfile', dumpfile ) % dump statistics for parsing
        cvx_begin SDP % Use CVX's SDP mode
        variable Z(dims, dims) symmetric
        variable t
        Z == semidefinite(dims);
        for i= 1:numRobots
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
        if(iter > 1)
            cost = cost + scale * trace(C*Z)
        end
        cost <= t
        minimize(t)
        cvx_end
        delete([dumpfile,'.mat'])
        cvx_solver_settings -clear dumpfile
        [V,D] = eig(Z);

        record = [record diag(D)];
        ratio = scale * trace(C*Z) / cost;
        if(ratio < tolerance)
            break;
        end
        U = V(:,1:end-3);
        C = U * U';
        scale = evaluate_objective_WYJ_by_X(Z, auxiliary_data) / (1 + trace(C*Z)) * 20;
    end
    mark_toc = toc(mark_tic);
    [U,S,V] = svd(Z);
    sqrt_S = sqrt(S);
    Yopt  = sqrt_S(1:3,1:3) * V(:,1:3)';
    if(det(Yopt(:,1:3)) < 0)
        Yopt = -Yopt;
    end
    R_opt_iteration = zeros(3, 3*numRobots);
    for k = 1:numRobots
        Y = Yopt(:,3*(k-1)+1:3*(k-1)+3);
        [U,S,V] = svd(Y);
        R_opt_iteration(:, 3*(k-1) + 1 : 3*(k-1) + 3) = U*V'*det(U*V');
    end
    f_X = evaluate_objective_WYJ_by_X(Z, auxiliary_data);
    g_RTR = evaluate_objective_WYJ(R_opt_iteration, auxiliary_data);

    result.C = C;
    result.scale = scale;
    result.R_opt = R_opt_iteration(:,1:3)'*R_opt_iteration;
    result.f_X = f_X;
    result.f_X_CX = cost;
    result.g_RTR = g_RTR;
    result.refuse_error = norm(Yopt'*Yopt - Z);
    result.rounding_error =  norm(R_opt_iteration - Yopt);
    result.l2error = norm(result.R_opt - data.data_gt.R_gt);
    result.time = mark_toc;
    result.record = record;
    result.X = Z;
    result.iteration = iter;
end