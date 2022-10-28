function egrad = Euclidean_gradient_WYJ(Y, problem_data)
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
            egrad_temp = B*X*A + A*X*B;
            egrad_temp = egrad_temp + egrad_temp';
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
                BX = e; BX(index_i,:) = X(index_j,:);
                XB = e; XB(:,index_j) = X(:,index_i);
%                 A = e; 
                BXA = e;
                AXB = e;
                for obsed = 1:numRobots
                    if(obsed == k)
                        continue;
                    end
                    I = information{k,obsed}(3*(i-1)+1:3*i,3*(j-1)+1:3*j);
%                     A(3*(obsed-1)+1:3*(obsed-1)+3, 3*(obsed-1)+1: 3*(obsed-1)+3) = ...
%                         information{k,obsed}(3*(i-1)+1:3*i,3*(j-1)+1:3*j);
                    BXA(index_i, 3*(obsed-1)+1:3*(obsed-1)+3) = BX(index_i, 3*(obsed-1)+1:3*(obsed-1)+3) * I; % 1 x 3
                    AXB(3*(obsed-1)+1:3*(obsed-1)+3, index_j) = I * XB(3*(obsed-1)+1:3*(obsed-1)+3, index_j); % 3 x 1
                end
%                 egrad_temp = B*X*A + A*X*B;
%                 egrad_temp = BX*A + A*XB;
                egrad_temp = AXB + BXA;
                egrad_temp = egrad_temp + egrad_temp';
                egrad = egrad + egrad_temp;
            end
        end
    end
end
egrad = egrad + problem_data.C + problem_data.C';
egrad = Y*egrad;
end

