function [cost, gradient] = evaluate_objective_WYJ(Y, problem_data)
% function [trQYtY, YQ] = evaluate_objective(Y, problem_data)
%
% This function computes and returns the value of the objective function
% tr(Ａ Y^T Y B Y^T Y) + tr(C Y^T Y ) + tr(D).  Optionally, it returns the product YQ as the second
% argument
numRobots =  problem_data.n;
X = Y'*Y;
if(~isfield(problem_data, 'information'))
    cost = vec(X)' * problem_data.BA * vec(X);
else
    e = zeros(3*numRobots,3*numRobots);
    cost = 0;
    information = problem_data.information;
    for k = 1:numRobots
        for i = 1:3
            for j = 1:3
                index_i = 3*(k-1)+i;
                index_j = 3*(k-1)+j;
                A = e; 
                M = X(:,index_i)*X(index_j,:);
                for obsed = 1:numRobots
                    if(obsed == k)
                        continue;
                    end
                    cost = cost + trace(information{k,obsed}(3*(i-1)+1:3*i,3*(j-1)+1:3*j)*M(3*(obsed-1)+1:3*(obsed-1)+3, 3*(obsed-1)+1: 3*(obsed-1)+3));
%                     cost = cost + X(index_j,3*(obsed-1)+1:3*(obsed-1)+3) * information{k,obsed}(3*(i-1)+1:3*i,3*(j-1)+1:3*j) * X(3*(obsed-1)+1:3*(obsed-1)+3,index_i);
                end
            end
        end
    end
end
cost = cost + trace(problem_data.C * X) + trace(problem_data.D);
end

