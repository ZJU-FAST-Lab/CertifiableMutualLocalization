function [log] = log_se3(T)
% @brief SE(3) matrix logarithm
%
% Equation is from Ethan Eade's reference: http://ethaneade.com/lie.pdf
% \f{align%}{
% \boldsymbol\omega &=\mathrm{skew\_offdiags}\Big(\frac{\theta}{2\sin\theta}(\mathbf R - \mathbf R^\top)\Big) \\
% \mathbf u &= \mathbf V^{-1}\mathbf t
% \f}
% where we have the following definitions
% \f{align%}{
% \theta &= \mathrm{arccos}((\mathrm{tr}(\mathbf R)-1)/2) \\
% \mathbf V^{-1} &= \mathbf I - \frac{1}{2} \lfloor \boldsymbol\omega \times\rfloor + \frac{1}{\theta^2}\Big(1-\frac{A}{2B}\Big)\lfloor
% \boldsymbol\omega \times\rfloor^2 \f}
%
% This function is based on the GTSAM one as the original implementation was a bit unstable.
% See the following:
% - https://github.com/borglab/gtsam/
% - https://github.com/borglab/gtsam/issues/746
% - https://github.com/borglab/gtsam/pull/780
%
w = log_so3(T(1:3, 1:3));
t = T(1:3, 4);
n = norm(w,2);
if (n < 1e-10)
    log = [w;t];
else
    W = skew(w / n);
    % Formula from Agrawal06iros, equation (14)
    % simplified with Mathematica, and multiplying in T to avoid matrix math
    Tan = tan(0.5 * n);
    Wt = W * t;
    u = t - (0.5 * n) * Wt + (1 - n / (2. * Tan)) * (W * Wt);
    log = [w;u];
end
end

