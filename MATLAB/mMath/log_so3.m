function [omega] = log_so3(R)
%SO(3) matrix logarithm
%
% This definition was taken from "Lie Groups for 2D and 3D Transformations" by Ethan Eade equation 17 & 18.
% http://ethaneade.com/lie.pdf
% \f{align%}{
% \theta &= \textrm{arccos}(0.5(\textrm{trace}(\mathbf{R})-1)) \\
% \lfloor\mathbf{v}\times\rfloor &= \frac{\theta}{2\sin{\theta}}(\mathbf{R}-\mathbf{R}^\top)
% @f}
%
% This function is based on the GTSAM one as the original implementation was a bit unstable.
% See the following:
% - https://github.com/borglab/gtsam/
% - https://github.com/borglab/gtsam/issues/746
% - https://github.com/borglab/gtsam/pull/780

% note switch to base 1
R11 = R(1, 1); R12 = R(1, 2); R13 = R(1, 3);
R21 = R(2, 1); R22 = R(2, 2); R23 = R(2, 3);
R31 = R(3, 1); R32 = R(3, 2); R33 = R(3, 3);
tr = trace(R);
%when trace == -1, i.e., when theta = +-pi, +-3pi, +-5pi, etc.
% we do something special
  if (tr + 1.0 < 1e-10) 
    if (abs(R33 + 1.0) > 1e-5)
      omega = (pi / sqrt(2.0 + 2.0 * R33)) * [R13; R23; 1.0 + R33];
    elseif (abs(R22 + 1.0) > 1e-5)
      omega = (pi / sqrt(2.0 + 2.0 * R22)) * [R12; 1.0 + R22; R32];
    else
      % if(std::abs(R.r1_.x()+1.0) > 1e-5)  This is implicit
      omega = (pi / sqrt(2.0 + 2.0 * R11)) * [1.0 + R11; R21; R31];
    end
  else
    tr_3 = tr - 3.0; % always negative
    if (tr_3 < -1e-7)
      theta = acos((tr - 1.0) / 2.0);
      magnitude = theta / (2.0 * sin(theta));
    else
      % when theta near 0, +-2pi, +-4pi, etc. (trace near 3.0)
      % use Taylor expansion: theta \approx 1/2-(t-3)/12 + O((t-3)^2)
      % see https:%github.com/borglab/gtsam/issues/746 for details
      magnitude = 0.5 - tr_3 / 12.0;
    end
    omega = magnitude * [R32 - R23; R13 - R31; R21 - R12];
   end
end

