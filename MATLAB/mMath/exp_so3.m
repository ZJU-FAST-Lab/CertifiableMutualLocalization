function [R] = exp_so3(w)
%SO(3) matrix exponential
%SO(3) matrix exponential mapping from the vector to SO(3) lie group.
% SO(3) matrix exponential mapping from the vector to SO(3) lie group.
% This formula ends up being the [Rodrigues formula](https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula).
% This definition was taken from "Lie Groups for 2D and 3D Transformations" by Ethan Eade equation 15.
% http://ethaneade.com/lie.pdf
% \f{align%}{
% \exp\colon\mathfrak{so}(3)&\to SO(3) \\
% \exp(\mathbf{v}) &=
% \mathbf{I}
% +\frac{\sin{\theta}}{\theta}\lfloor\mathbf{v}\times\rfloor
% +\frac{1-\cos{\theta}}{\theta^2}\lfloor\mathbf{v}\times\rfloor^2 \\
% \mathrm{where}&\quad \theta^2 = \mathbf{v}^\top\mathbf{v}
 % @f}
 w_x = skew(w);
 theta = norm(w,2);
 if(theta < 1e-7)
     A = 1; 
     B = 0.5;
 else
     A = sin(theta) / theta;
     B = (1 - cos(theta)) / (theta * theta);
 end
 if (theta == 0)
    R = eye(3);
      else
        R = eye(3) + A * w_x + B * w_x * w_x;
 end
end

