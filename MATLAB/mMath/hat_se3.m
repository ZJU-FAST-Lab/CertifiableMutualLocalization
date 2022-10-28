function [mat] = hat_se3(vec)
% @brief Hat operator for R^6 -> Lie Algebra se(3)
%
% \f{align%}{
% \boldsymbol\Omega^{\wedge} = \begin{bmatrix} \lfloor \boldsymbol\omega \times\rfloor & \mathbf u \\ \mathbf 0 & 0 \end{bmatrix}
% \f}

mat = zeros(4,4);
mat(1:3,1:3) = skew(vec(1:3));
mat(1:3,4) = vec(4:6);

end

