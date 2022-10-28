function w = vee(w_x)
% @brief Returns vector portion of skew-symmetric
%
% See skew() for details.

w = [w_x(3,2); w_x(1,3); w_x(2,1)];

end
