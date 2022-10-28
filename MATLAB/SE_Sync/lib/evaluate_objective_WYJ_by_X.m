function cost = evaluate_objective_WYJ_by_X(X, problem_data)
% function [trQYtY, YQ] = evaluate_objective(Y, problem_data)
%
% This function computes and returns the value of the objective function
% tr(Ａ Y^T Y B Y^T Y) + tr(C Y^T Y ) + tr(D).  Optionally, it returns the product YQ as the second
% argument


numRobots =  problem_data.n;
if(~isfield(problem_data, 'information'))
    cost = vec(X)' * problem_data.BA * vec(X);

else
    information = problem_data.information;
    Q = zeros(9*numRobots*numRobots, 9*numRobots*numRobots);
    for k = 1:numRobots
        indexk_min = 9*numRobots * (k-1) + 1;
        indexk_max = 9*numRobots * k;
        G_k = zeros(9*numRobots, 9*numRobots);
        for i = 1:3
            for j = 1:3
%                 index_i = 3*(k-1)+i;
%                 index_j = 3*(k-1)+j;
                G_k_obsed = zeros(3*numRobots, 3*numRobots);
                for obsed = 1:numRobots
                    if(obsed == k)
                        continue;
                    end
                    index_obsed_min = 3*(obsed-1) + 1;
                    index_obsed_max = 3*(obsed-1) + 3;
                    G_k_obsed(index_obsed_min:index_obsed_max, index_obsed_min:index_obsed_max) = information{k,obsed}(3*(i-1)+1:3*i,3*(j-1)+1:3*j);
                end
                index_i_min = 3*numRobots * (i-1) + 1;
                index_i_max = 3*numRobots * i;
                index_j_min = 3*numRobots * (j-1) + 1;
                index_j_max = 3*numRobots * j;
                G_k(index_i_min:index_i_max, index_j_min:index_j_max) = G_k_obsed;
            end
        end
        Q(indexk_min:indexk_max, indexk_min:indexk_max) = G_k;
    end
    cost = vec(X)' * Q * vec(X);

%     e = zeros(3*numRobots,3*numRobots);
%     cost = 0;
%     information = problem_data.information;
%     for k = 1:numRobots
%         for i = 1:3
%             for j = 1:3
%                 index_i = 3*(k-1)+i;
%                 index_j = 3*(k-1)+j;
%                 M = X(:,index_i)*X(index_j,:);
%                 for obsed = 1:numRobots
%                     if(obsed == k)
%                         continue;
%                     end
%                     cost = cost + trace(information{k,obsed}(3*(i-1)+1:3*i,3*(j-1)+1:3*j)*M(3*(obsed-1)+1:3*(obsed-1)+3, 3*(obsed-1)+1: 3*(obsed-1)+3));
%                 end
%             end
%         end
%     end
end
cost = cost + trace(problem_data.C * X) + trace(problem_data.D);

end

