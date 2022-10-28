function [R,f] = project_rotation( M, fmin )
% R = project_rotation( M )
% 
% Takes a general linear matrix M and projects it to
% the rotation matrix in SO(3) minimizing the chordal distance.
% This is equivalent to find the closest matrix under Frobenius norm
% (see "Hartley et al. Rotation averaging" in IJCV).
% 
% The problem is solved using the SVD decomposition of M:
%   M = U*S*V'
% If det(U*V')>0, R = U*V'
% BUT if det(U*V')<0, R = U*V'
% 
% 
% [R,f] = project_rotation( M, fmin )
% If the projected rotation seeks to minimize a certain objective fmin,
% it may be preferrable to search for the rotation
%   R=U*S*V', S=diag(s_i), s_i=+-1
% that minimizes fmin (this can be done by brute-force search)

% Project to valid rotation
[U,D,V] = svd(M);
if nargin == 1
  if det(U*V')>0
    R = U*V';
  else
    R = U*diag([1 1 -1])*V';
  end
else
  if det(U*V')>0
    signs = {[1 1 1],[1 -1 -1],[-1 1 -1],[-1 -1 1]};
  else
    signs = {[-1 1 1],[1 -1 1],[1 1 -1],[-1 -1 -1]};
  end
  numCombs = numel(signs);
  fvalues = zeros(1,numCombs);
  for i=1:numCombs
    fvalues(i) = fmin(U*diag(signs{i})*V');
  end
  [f,bestIdx] = min(fvalues);
  R = U*diag(signs{bestIdx})*V';
end
    
