function Lambda_blocks = compute_Lambda(Yopt, problem_data)
%function Lambda_blocks = compute_Lambda(Yopt, problem_data)
%
% Given an estimated local minimum Yopt for the (possibly lifted) relaxation,
% this function computes and returns the block-diagonal elements of the
% corresponding Lagrange multiplier 

Lambda_blocks = zeros(problem_data.d, problem_data.d*problem_data.n);

QYt = Qproduct(Yopt', problem_data);

for k = 1:problem_data.n
    imin = problem_data.d*(k-1) + 1;
    imax = problem_data.d*(k-1) + problem_data.d;
    
    B = QYt(imin:imax, :) * Yopt(:, imin:imax);
    Lambda_blocks(:, imin:imax) = .5 * (B + B');
end

