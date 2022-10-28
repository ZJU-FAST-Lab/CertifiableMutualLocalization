function [result] = local_optimization(data)

    auxiliary_data = data.auxiliary_data;
    numRobots = data.numRobots;
 
    %% Rieman Optimization Refine
    %% Set Manopt options
    manopt_options.tolgradnorm = 1e-2;  % Stopping tolerance for norm of Riemannian gradient
    manopt_options.miniter = 1;  % Minimum number of outer iterations (i.e. accepted update steps) to perform
    manopt_options.maxiter = 100;  % Maximum number of outer iterations (i.e. accepted update steps) to perform
    manopt_options.maxinner = 500;  % Maximum number of iterations for the conjugate-gradient method used to compute approximate Newton steps
    
    %% Set SE-Sync options
    se_sync_opts.r0 = 3;  % Initial maximum-rank parameter at which to start the Riemannian Staircase
    se_sync_opts.rmax = se_sync_opts.r0;  % Maximum maximum-rank parameter at which to terminate the Riemannian Staircase
    se_sync_opts.eig_comp_rel_tol = 1e-5;  % Relative tolerance for the minimum-eigenvalue computation used to test for second-order optimality with MATLAB's eigs() function
    se_sync_opts.min_eig_threshold = -1e-2;  % Minimum eigenvalue threshold for accepting a maxtrix as numerically positive-semidefinite
    se_sync_opts.relative_func_decrease_tol = 1e-6;  % Additional stopping criterion for Manopt: stop if the relative function decrease between two successive accepted iterates is less than this value


    %% Run SE-Sync
    t_start = tic();
    [SDPval, Yopt, xhat, Fxhat, se_sync_info, auxiliary_data] = SE_Sync_WYJ(auxiliary_data, manopt_options, se_sync_opts);
    t_end = toc(t_start);

    if(det(Yopt(:,1:3)) < 0)
        Yopt = -Yopt;
    end
    R_opt_manifold = zeros(3, 3*numRobots);
    for k = 1:numRobots
        Y = Yopt(:,3*(k-1)+1:3*(k-1)+3);
        [U,S,V] = svd(Y);
        R_opt_manifold(:, 3*(k-1) + 1 : 3*(k-1) + 3) = U*V'*det(U*V');
    end
    g_RTR = evaluate_objective_WYJ(R_opt_manifold, auxiliary_data);
    f_X = evaluate_objective_WYJ_by_X(Yopt'*Yopt, auxiliary_data);

    result.R_opt = R_opt_manifold(:,1:3)' * R_opt_manifold;
    result.f_X = f_X;
    result.g_RTR = g_RTR;
    result.refuse_error = 0.0;
    result.rounding_error = norm(R_opt_manifold - Yopt);
    result.time = t_end;
    result.l2error = norm(result.R_opt - data.data_gt.R_gt);
    result.X = Yopt'*Yopt;
end