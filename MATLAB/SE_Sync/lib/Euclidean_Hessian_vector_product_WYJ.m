function Hvec = Euclidean_Hessian_vector_product_WYJ(Y, Ydot, problem_data)
%function Hvec = Euclidean_Hessian_vector_product(Y, Ydot, problem_data)
%
% This function computes and returns the value of the Euclidean Hessian at
% the point Y evaluated along the tangent direction Ydot.

numRobots =  problem_data.n;
X = Y'*Y;
X1 = Y'*Ydot;
e = zeros(3*numRobots,3*numRobots);
Hvec1 = zeros(3*numRobots,3*numRobots);
Hvec2 = zeros(3*numRobots,3*numRobots);

if(~isfield(problem_data, 'information'))
    Q = problem_data.BA;
    for i = 1:3*numRobots
        for j = i:3*numRobots
            B = e; B(i,j)=1;
            A = Q((j-1)*3*numRobots+1:j*3*numRobots, (i-1)*3*numRobots+1: i*3*numRobots);
            Hvec_temp1 = B*X1*A + B'*X1*A'+ A'*X1*B' + A*X1*B;
            Hvec_temp1 = Hvec_temp1 + Hvec_temp1';
            Hvec_temp2 = B*X*A + B'*X*A'+ A'*X*B' + A*X*B;
            if(i ~= j)
                Hvec_temp1 = Hvec_temp1 * 2;
                Hvec_temp2 = Hvec_temp2 * 2;
            end
            Hvec1 = Hvec1 + Hvec_temp1;
            Hvec2 = Hvec2 + Hvec_temp2;
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
    
                BX1 = e; BX1(index_i,:) = X1(index_j,:);
                X1B = e; X1B(:,index_j) = X1(:,index_i);
                BX = e; BX(index_i,:) = X(index_j,:);
                XB = e; XB(:,index_j) = X(:,index_i);

                BXA = e;
                AXB = e;
                BX1A = e;
                AX1B = e;
                for obsed = 1:numRobots
                    if(obsed == k)
                        continue;
                    end
                    I = information{k,obsed}(3*(i-1)+1:3*i,3*(j-1)+1:3*j);
%                     A(3*(obsed-1)+1:3*(obsed-1)+3, 3*(obsed-1)+1: 3*(obsed-1)+3) = ...
%                         information{k,obsed}(3*(i-1)+1:3*i,3*(j-1)+1:3*j);

                    BXA(index_i, 3*(obsed-1)+1:3*(obsed-1)+3) = BX(index_i, 3*(obsed-1)+1:3*(obsed-1)+3) * I; % 1 x 3
                    AXB(3*(obsed-1)+1:3*(obsed-1)+3, index_j) = I * XB(3*(obsed-1)+1:3*(obsed-1)+3, index_j); % 3 x 1
                    BX1A(index_i, 3*(obsed-1)+1:3*(obsed-1)+3) = BX1(index_i, 3*(obsed-1)+1:3*(obsed-1)+3) * I; % 1 x 3
                    AX1B(3*(obsed-1)+1:3*(obsed-1)+3, index_j) = I * X1B(3*(obsed-1)+1:3*(obsed-1)+3, index_j); % 3 x 1
                end
        
%                 Hvec_temp1 = BX1*A  + A*X1B;
%                 Hvec_temp2 = BX*A +  A*XB;
                Hvec_temp1 = BX1A  + AX1B;
                Hvec_temp2 = BXA +  AXB;
                Hvec_temp1 = 2*(Hvec_temp1 + Hvec_temp1');
                Hvec_temp2 = Hvec_temp2 + Hvec_temp2';
    
                Hvec1 = Hvec1 + Hvec_temp1;
                Hvec2 = Hvec2 + Hvec_temp2;
               
            end
        end
    end
end
Hvec = Y*Hvec1+Ydot*Hvec2;
Hvec = Hvec + Ydot*(problem_data.C + problem_data.C');
end

