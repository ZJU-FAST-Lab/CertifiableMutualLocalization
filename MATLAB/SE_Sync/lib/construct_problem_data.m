function problem_data = construct_problem_data(measurements)
%function problem_data = construct_problem_data(measurements)
%
% This helper function accepts a MATLAB struct containing the raw 
% measurements specifying a special Euclidean synchronization problem, and 
% returns another struct containing the data matrices required by the 
% SE-Sync algorithm.  Formally:
%
% INPUT: A MATLAB struct 'measurements' containing the following fields
% (see eq. (11) in the long-form version of the paper for details):
% edges:  An (mx2)-dimension encoding the edges in the measurement network;
%     edges(k, :) = [i,j] means that the kth measurement is of the
%     relative transform from pose i to pose j.  NB:  This indexing scheme 
%     requires that the states x_i are numbered sequentially as 
%     x_1, ... x_n.
% R:  An m-dimensional cell array whose kth element is the rotational part
%     of the kth measurement
% t:  An m-dimensional cell array whose kth element is the translational
%     part of the kth measurement
% kappa:  An m-dimensional cell array whose kth element gives the precision
%     of the rotational part of the kth measurement. 
% tau:  An m-dimensional cell array whose kth element gives the precision
%     of the translational part of the kth measurement.
%
% 
%
% OUTPUT:  A MATLAB struct 'problem_data' containing various data matrices
% that the SE-Sync algorithm requires:
%
% d:  The dimensional parameter of the special Euclidean group over which 
%     the estimation takes place (typically d = 2 or 3).
% n:  The number of states to be estimated.
% m:  The number of available measurements.
% ConLap:  The connection Laplacian for the rotational measurements (see
%     eq. (15) in the paper.
% A:  The oriented incidence matrix for the underlying directed graph of
%     measurements (see eq. (7) in the paper).
% Ared:  The reduced incidence matrix obtained by removing the final 
%     row of A.
% L:  The lower-triangular factor of a thin LQ decomposition of the reduced
%     incidence matrix Ared (see eq. (40) in the paper).
% Omega:  The diagonal matrix of translational matrix precisions (see eq.
%     (23) in the paper).
% T:  The sparse matrix of translational observations definedin equation 
%     (24) in the paper.
% V:  The matrix of translational observations defined in equation 
%     (16) in the paper

%Given the 'measurements' struct returned by 'readG2oDataset3D', this
%function constructs and returns the data matrices defining the pose-graph
%relaxation problem

% Set additional variables
problem_data.d = length(measurements.t{1});
problem_data.n = max(max(measurements.edges));
problem_data.m = size(measurements.edges, 1);


% Construct connection Laplacian for the rotational measurements
tic();
ConLap = construct_connection_Laplacian(measurements);
problem_data.ConLap = ConLap;
t = toc();
disp(sprintf('Constructed rotational connection Laplacian in %g seconds', t));

% Construct the oriented incidence matrix for the underlying directed graph
% of measurements
tic();
A = construct_incidence_matrix(measurements);
problem_data.A = A;
t = toc();
disp(sprintf('Constructed oriented incidence matrix in %g seconds', t));

% Construct the reduced oriented incidence matrix
Ared = A(1:problem_data.n-1, :);
problem_data.Ared = Ared;

% Construct the lower-triangular factor for the reduced oriented incidence
% matrix
tic();
reducedLaplacian = Ared * Ared';
L = chol(reducedLaplacian, 'lower');
problem_data.L = L;
t = toc();
disp(sprintf('Computed lower-triangular factor of reduced incidence matrix in %g seconds', t));

tic();
[T, Omega] = construct_translational_matrices(measurements);
V = construct_V_matrix(measurements);
t = toc();
disp(sprintf('Constructed translational observation and measurement precision matrices in %g seconds', t));

problem_data.T = T;
problem_data.Omega = Omega;
problem_data.V = V;
% B = Ared' * inv(L') * inv(L) * Ared;
% I = eye(size(B));
% Pi = I - B;
% problem_data.Q = ConLap + T' * Pi * Omega * Pi * T;

end

