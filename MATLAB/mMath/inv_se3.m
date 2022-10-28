function [inv_T] = inv_se3(T)
inv_T = eye(4);
inv_T(1:3,1:3) = T(1:3,1:3)';
inv_T(1:3,4) = -inv_T(1:3,1:3) * T(1:3,4);
end

