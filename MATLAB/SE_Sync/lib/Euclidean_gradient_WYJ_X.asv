function egrad = Euclidean_gradient_WYJ_X(Y, problem_data)
%
% This function computes and returns the value of the Euclidean gradient of
% the objective function: nabla F(Y) = 2YQ.


numRobots =  problem_data.n;
X = Y'*Y;
e = zeros(3*numRobots,3*numRobots);
egrad = zeros(3*numRobots,3*numRobots);
if(~isfield(problem_data, 'information'))
    Q = problem_data.BA;
    for i = 1:3*numRobots
        for j = i:3*numRobots
            B = e; B(i,j)=1;
            A = Q((j-1)*3*numRobots+1:j*3*numRobots, (i-1)*3*numRobots+1: i*3*numRobots);
            egrad_temp = A'*X*B' + A*X*B;
            if(i ~= j)
                egrad_temp = egrad_temp * 2;
            end
            egrad = egrad + egrad_temp;
        end
    end
else
    information = problem_data.information;
    for k = 1:numRobots
        for i = 1:3
            for j = 1:3
                index_i = 3*(k-1)+i;
                index_j = 3*(k-1)+j;
                B = e; B(index_i, index_j) = 1;
                A = e; 
                for obsed = 1:numRobots
                    if(obsed == k)
                        continue;
                    end
                    A(3*(obsed-1)+1:3*(obsed-1)+3, 3*(obsed-1)+1: 3*(obsed-1)+3) = ...
                        information{k,obsed}(3*(i-1)+1:3*i,3*(j-1)+1:3*j);
                end
                egrad_temp = A'*X*B' + A*X*B;
                egrad = egrad + egrad_temp;
            end
        end
    end
end
egrad = egrad + problem_data.C';
end

