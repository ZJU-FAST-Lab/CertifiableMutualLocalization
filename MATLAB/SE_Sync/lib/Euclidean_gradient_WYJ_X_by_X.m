function egrad = Euclidean_gradient_WYJ_X_by_X(X, problem_data)
%
% This function computes and returns the value of the Euclidean gradient of
% the objective function: nabla F(Y) = 2YQ.


numRobots =  problem_data.n;
e = zeros(3*numRobots,3*numRobots);
egrad = zeros(3*numRobots,3*numRobots);
Q = problem_data.BA;
for i = 1:3*numRobots
    for j = 1:3*numRobots
        B = e; B(i,j)=1;
        A = Q((j-1)*3*numRobots+1:j*3*numRobots, (i-1)*3*numRobots+1: i*3*numRobots);
        egrad_temp = A'*X*B' + A*X*B;
        egrad = egrad + egrad_temp;
    end
end
egrad = egrad + problem_data.C';
end

