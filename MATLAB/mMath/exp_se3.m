function [T] = exp_se3(vec)
 % @brief SE(3) matrix exponential function
 %
 % Equation is from Ethan Eade's reference: http://ethaneade.com/lie.pdf
 % \f{align%}{
 % \exp([\boldsymbol\omega,\mathbf u])&=\begin{bmatrix} \mathbf R & \mathbf V \mathbf u \\ \mathbf 0 & 1 \end{bmatrix} \\[1em]
 % \mathbf R &= \mathbf I + A \lfloor \boldsymbol\omega \times\rfloor + B \lfloor \boldsymbol\omega \times\rfloor^2 \\
 % \mathbf V &= \mathbf I + B \lfloor \boldsymbol\omega \times\rfloor + C \lfloor \boldsymbol\omega \times\rfloor^2
 % \f}
 % where we have the following definitions
 % \f{align%}{
 % \theta &= \sqrt{\boldsymbol\omega^\top\boldsymbol\omega} \\
 % A &= \sin\theta/\theta \\
 % B &= (1-\cos\theta)/\theta^2 \\
 % C &= (1-A)/\theta^2
 % \f}

% Precompute our values
  w = vec(1:3);
  u = vec(4:6);
  theta = norm(w,2);
  wskew = skew(w);

  % Handle small angle values
  if (theta < 1e-7)
    A = 1;
    B = 0.5;
    C = 1.0 / 6.0;
  else
    A = sin(theta) / theta;
    B = (1 - cos(theta)) / (theta * theta);
    C = (1 - A) / (theta * theta);
  end

  % Matrices we need V and Identity
  I_33 = eye(3);
  V = I_33 + B * wskew + C * wskew * wskew;

  % Get the final matrix to return
  T = zeros(4);
  T(1:3,1:3) = I_33 + A * wskew + B * wskew * wskew;
  T(1:3,4) = V * u;
  T(4, 4) = 1;
end

