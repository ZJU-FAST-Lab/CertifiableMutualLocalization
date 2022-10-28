function M = offsym(v)
% M = offsym(v)
% Convert 3-element list into its offdiagonal-symmetric counterpart
% This is an equivalent concept to skew matrix [-]_x of vectors
% 
% The offsym operator is commutative with right vector product:
%   offsym(a)*b = offsym(b)*a

% This only makes sense for 3-dimensional vectors or lists
assert(numel(v)==3)

if ~iscell(v)
  M = [ 0   v(3)  v(2)
       v(3)  0    v(1)
       v(2)  v(1)  0  ];
else
  s  = size(v{1});
  Os = zeros(s);
  M  = [  Os   v{3} v{2};
         v{3}   Os  v{1};
         v{2}  v{1}  Os ];
end
end