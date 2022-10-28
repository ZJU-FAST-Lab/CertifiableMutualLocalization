function [SDPval, Yopt, xhat, Fxhat, SE_Sync_info, auxiliary_data] = SE_Sync_Q(auxiliary_data, manopt_options, SE_Sync_opts, Y0)
%function [SDPval, Yopt, xhat, Fxhat, SE_Sync_info, auxiliary_data] = SE_Sync(measurements, manopt_options, SE_Sync_opts, Y0)
%
% This function implements the SE-Sync algorithm described in the paper
% "SE-Sync: A Certifiably Correct Algorithm for Synchronization over the
% Special Euclidean Group"
% 
% INPUTs:
% 
% measurements:  A MATLAB struct containing the data describing the special
%   Euclidean synchronization problem (see eq. (11) in the paper for 
%   details). Specifically, measurements must contain the following fields:
%   edges:  An (mx2)-dimension encoding the edges in the measurement 
%     network; edges(k, :) = [i,j] means that the kth measurement is of the
%     relative transform x_i^{-1} x_j.  NB:  This indexing scheme requires
%     that the states x_i are numbered sequentially as x_1, ... x_n.
%   R:  An m-dimensional cell array whose kth element is the rotational part
%     of the kth measurement
%   t:  An m-dimensional cell array whose kth element is the translational
%     part of the kth measurement
%   kappa:  An m-dimensional cell array whose kth element gives the 
%     precision of the rotational part of the kth measurement. 
%   tau:  An m-dimensional cell array whose kth element gives the precision
%     of the translational part of the kth measurement.
%
% manopt_options:  A MATLAB struct containing various options that 
%       determine the behavior of Manopt's Riemannian truncated-Newton 
%       trust-region method, which we use to solve instances of the 
%       rank-restricted form of the semidefinite relaxation.  This struct 
%       contains the following fields (among others, see the Manopt
%       documentation)
%   tolgradnorm:  Stopping criterion; norm tolerance for the Riemannian gradient
%   maxinner:  Maximum number of Hessian-vector products to evaluate as part
%      of the truncated conjugate-gradient procedure used to compute update
%      steps.
%   miniter:  Minimum number of outer iterations (update steps).
%   maxiter:  Maximum number of outer iterations (update steps).
%   maxtime:  Maximum permissible elapsed computation time (in seconds).
%
% SE_Sync_opts:  A MATLB struct determining the behavior of the SE-Sync
%       algorithm.  This struct contains the following fields:
%   r0:  The initial value of the maximum-rank parameter r at which to 
%      start the Riemannian Staircase
%   rmax:  The maximum value of the maximum-rank parameter r.
%   eig_comp_rel_tol:  Relative tolerance for the minimum-eigenvalue 
%      computation needed to verify second-order optimality using MATLAB's 
%      eigs command (typical values here are on the order of 10^-5)
%   min_eig_threshold:  Threshold for the minimum eigenvalue in order to
%      consider the matrix Q - Lambda to be positive semidefinite.  Typical
%      values here should be small-magnitude negative numbers, e.g. 10^-4
%   rel_func_decrease_tol:  An additional stopping criterion for the Manopt
%     solver.  Terminate whenever the relative decrease in function value
%     between subsequenct iterations is less than this value (in the range
%     (0,1) ).
%
% Y0:  [Optional]  An initial point on the manifold St(d, r)^n at which to
%      initialize the first Riemannian optimization problem.  If this
%      parameter is not passed, a randomly-sampled point is used instead.
%
%
% OUTPUTS:
%
% SDPval:  The optimal value of the semidefinite relaxation
% Yopt:  A symmetric factor of an optimal solution Zopt = Yopt' * Yopt for
%      the semidefinite relaxation.
% xhat: A struct containing the estimate for the special Euclidean
%   synchronization problem.  It has the following two fields:
%   Rhat:  A d x dn matrix whose (dxd)-block elements give the rotational 
%   state estimates.
%   that: a d x n matrix whose columsn give the translational state estimates.
% Fxhat:  The objective value of the rounded solution xhat.
%
% SE_Sync_info:  A MATLAB struct containing various possibly-interesting
%   bits of information about the execution of the SE-Sync algorithm.  The
%   fields are:
%   optimization_times:  A vector containing the elapsed computation times 
%     for solving the optimization problem at each level of the Riemannian
%     Staircase.
%   SDPLRvals:  A vector containing the optimal value of the optimization
%     problem solved at each level of the Riemannian Staircase
%   min_eig_times:  A vector containing the elapsed computation times for
%     performing the minimum-eigenvalue computation necessary to check for
%     optimality of Yopt as a solution of the SDP after solving each
%     Riemannian optimization problem to first-order.
%   min_eig_vals:  A vector containing the corresponding minimum
%      eigenvalues.
%   total_computation_time:  The elapsed computation time of the complete
%      SE-Sync algorithm
%   manopt_info:  The info struct returned by the Manopt solver for the
%      during its last execution (i.e. when solving the last explored level
%      the Riemannian Staircase).
%
% auxiliary_data:  A MATLAB struct containing several auxiliary matrices
% constructed from the input measurements that are used internally
% throughout the SE-Sync algorithm.  Specifically, this struct contains the
% following fields:
%   n:  The number of group elements (poses) to estimate
%   m:  The number of relative measurements
%   d:  The dimension of the Euclidean space on which these group elements
%       act (generally d is 2 or 3).
%   ConLap:  The connection Laplacian for the set of rotational
%       measurements; see eq. (15) in the paper.
%   A:  An oriented incidence matrix for the directed graph of
%       measurements; see eq. (7) in the paper
%   Ared:  The reduced oriented incidence matrix obtained by removing the
%       final row of A.
%   L:  A sparse lower-triangular factor of a thin LQ decomposition of 
%       Ared; see eq. (40) in the paper
%   T:  The sparse matrix of translational observations defined in eq. (24)
%       in the paper
%   Omega:  The diagonal matrix of translational measurement precisions;
%       defined in eq. (23).
%   V:  The sparse translational data matrix defined in eq. (16) in the
%       paper.


disp(sprintf('\n\n========== SE-Sync ==========\n'));

timerVal = tic();


%% Construct problem data matrices from input
disp(sprintf('INITIALIZATION:\n'));
 
%% INITIALIZATION 

% The maximum number of levels in the Riemannian Staircase that we will
% need to explore
max_num_iters = SE_Sync_opts.rmax - SE_Sync_opts.r0 + 1;

% Allocate storage for state traces
optimization_times = zeros(1, max_num_iters);
SDPLRvals = zeros(1, max_num_iters);
min_eig_times = zeros(1, max_num_iters);
min_eig_vals = zeros(1, max_num_iters);


% Set up Manopt problem
% Check if a solver was explicitly supplied
if(~isfield(manopt_options, 'solver'))
    % Use the trust-region solver by default
    manopt_options.solver = @trustregions;
end
solver_name = func2str(manopt_options.solver);
if (~strcmp(solver_name, 'trustregions') && ~strcmp(solver_name, 'conjugategradient') && ~strcmp(solver_name, 'steepestdescent'))
    error(sprintf('Unrecognized Manopt solver: %s', solver_name));
end
disp(sprintf('Solving Riemannian optimization problems using Manopt''s "%s" solver\n', solver_name));

% Set cost function handles
manopt_data.cost = @(Y) evaluate_objective(Y', auxiliary_data);
manopt_data.egrad = @(Y) Euclidean_gradient(Y', auxiliary_data)';
manopt_data.ehess = @(Y, Ydot) Euclidean_Hessian_vector_product(Y', Ydot', auxiliary_data)';

% We optimize over the manifold M := St(d, r)^N, the N-fold product of the 
% (Stiefel) manifold of orthonormal d-frames in R^r.
manopt_data.M = stiefelstackedfactory(auxiliary_data.n, auxiliary_data.d, SE_Sync_opts.r0);

% Set additional stopping criterion for Manopt: stop if the relative
% decrease in function value between successive iterates drops below the
% threshold specified in SE_Sync_opts.relative_func_decrease_tol
if(strcmp(solver_name, 'trustregions'))
    manopt_options.stopfun = @(manopt_problem, x, info, last) relative_func_decrease_stopfun(manopt_problem, x, info, last, SE_Sync_opts.relative_func_decrease_tol);
end



% Check if an initial point was supplied
if nargin < 4
    disp(sprintf('Initializing Riemannian Staircase with randomly-sampled initial point on St(%d,%d)^%d\n', auxiliary_data.d, SE_Sync_opts.r0, auxiliary_data.n));
    % Sample a random point on the Stiefel manifold as an initial guess
    Y0 = manopt_data.M.rand()';
else
    disp(sprintf('Using user-supplied initial point Y0 in Riemannian Staircase\n'));
end

% Counter to keep track of how many iterations of the Riemannian Staircase
% have been performed
iter = 0;   

%%  RIEMANNIAN STAIRCASE
for r = SE_Sync_opts.r0 : SE_Sync_opts.rmax
    iter = iter + 1;  % Increment iteration number
    
    % Starting at Y0, use Manopt's truncated-Newton trust-region method to
    % descend to a first-order critical point.
    
    disp(sprintf('RIEMANNIAN STAIRCASE (level r = %d):', r)); 
    
    [YoptT, Fval, manopt_info, manopt_options] = manoptsolve(manopt_data, Y0', manopt_options);
    Yopt = YoptT';
    SDPLRval = Fval(end);
    
    % Store the optimal value and the elapsed computation time
    SDPLRvals(iter) = SDPLRval;
    optimization_times(iter) = manopt_info(end).time;

    
    % Augment Yopt by padding with an additional row of zeros; this
    % preserves Yopt's first-order criticality while ensuring that it is
    % rank-deficient
    
    Yplus = vertcat(Yopt, zeros(1, auxiliary_data.d * auxiliary_data.n));
       
    
    disp(sprintf('\nChecking second-order optimality...'));
    % At this point, Yplus is a rank-deficient critial point, so check
    % 2nd-order optimality conditions
    
    % Compute Lagrange multiplier matrix Lambda corresponding to Yplus
    Lambda = compute_Lambda(Yopt, auxiliary_data);
    
    % Compute minimum eigenvalue/eigenvector pair for Q - Lambda
    tic();
    [lambda_min, v] = Q_minus_Lambda_min_eig(Lambda, auxiliary_data, SE_Sync_opts.eig_comp_rel_tol);
    min_eig_comp_time = toc();
    
    % Store the minimum eigenvalue and elapsed computation times
    min_eig_vals(iter) = lambda_min;
    min_eig_times(iter) = min_eig_comp_time;
    
    if( lambda_min > SE_Sync_opts.min_eig_threshold)
        % Yopt is a second-order critical point
        disp(sprintf('Found second-order critical point! (minimum eigenvalue = %g, elapsed computation time %g seconds)\n', lambda_min, min_eig_comp_time));
        break;
    else
        disp(sprintf('Saddle point detected (minimum eigenvalue = %g,  elapsed computation time %g seconds)', lambda_min, min_eig_comp_time));
        % lambda_min is a negative eigenvalue of Q - Lambda, so the KKT
        % conditions for the semidefinite relaxation are not satisfied;
        % this implies that Yplus is a saddle point of the rank-restricted
        % semidefinite optimization.  Fortunately, the eigenvector v
        % corresponding to lambda_min can be used to provide a descent
        % direction from this saddle point, as described in Theorem 3.9 of
        % the paper "A Riemannian Low-Rank Method for Optimization over
        % Semidefinite Matrices with Block-Diagonal Constraints".
        
        % Define the vector Ydot := e_{r+1} * v'; this is tangent to the
        % manifold St(d, r+1)^n at Yplus and provides a direction of
        % negative curvature
        disp('Computing escape direction...');
        Ydot = vertcat(zeros(r, auxiliary_data.d * auxiliary_data.n), v');
        
        % Compute the directional derivative of F at Yplus along Ydot
        dF0 = trace(Euclidean_gradient(Yplus, auxiliary_data)*Ydot');
        if dF0 > 0
            Ydot = -Ydot;
        end
        
        % Augment the dimensionality of the Stiefel manifolds in
        % preparation for the next iteration
        
        manopt_data.M = stiefelstackedfactory(auxiliary_data.n, auxiliary_data.d, r+1);
        
        % Perform line search along the escape direction Ydot to escape the
        % saddle point and obtain the initial iterate for the next level in
        % the Staircase
        
        disp('Line searching along escape direction to escape saddle point...');
        tic();
%         [stepsize, Y0T] = linesearch_decrease(manopt_data, Yplus', Ydot', SDPLRval);
        [stepsize, Y0T] = linesearch_vector(manopt_data, Yplus', Ydot', SDPLRval);
        line_search_time = toc();
        Y0 = Y0T';
        disp(sprintf('Line search completed (elapsed computation time %g seconds)\n', line_search_time));
    end
end

disp(sprintf('===== END RIEMANNIAN STAIRCASE =====\n'));

%% POST-PROCESSING

% Return optimal value of the SDP (in the case that a rank-deficient,
% second-order critical point is obtained, this is equal to the optimum
% value obtained from the Riemannian optimization

SDPval = SDPLRval;

disp('Rounding solution...');
% Round the solution
tic();
[Rhat, singular_values, determinants] = round_solution(Yopt, auxiliary_data);
solution_rounding_time = toc();
disp(sprintf('Elapsed computation time: %g seconds\n', solution_rounding_time));

disp('Recovering translational estimates...');
% Recover the optimal translational estimates
tic();
% that = recover_translations(Rhat, auxiliary_data);
translation_recovery_time = toc();
disp(sprintf('Elapsed computation time: %g seconds\n', translation_recovery_time));

xhat.R = Rhat;
% xhat.t = that;

Fxhat = evaluate_objective(Rhat, auxiliary_data);

disp(sprintf('Suboptimality bound of recovered solution xhat: %g\n', Fxhat - SDPval));
total_computation_time = toc(timerVal);

disp(sprintf('Total elapsed computation time: %g seconds\n', total_computation_time));

% Output info
SE_Sync_info.SDPLRvals = SDPLRvals(1:iter);
SE_Sync_info.optimization_times = optimization_times(1:iter);
SE_Sync_info.min_eig_vals = min_eig_vals(1:iter);
SE_Sync_info.min_eig_times = min_eig_times(1:iter);
SE_Sync_info.manopt_info = manopt_info;
SE_Sync_info.total_computation_time = total_computation_time;
SE_Sync_info.determinants = determinants;

disp(sprintf('\n===== END SE-SYNC =====\n'));

        
        
        
        
        

    
    
    
    

    
    



end

