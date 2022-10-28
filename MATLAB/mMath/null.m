function Z = null(A,tol)
% Z = null(A)
% Z = null(A,tol)
% 
% Z = NULL(A) is an orthonormal basis for the null space of A obtained
% from the singular value decomposition.
% Namely, all those right singular vectors with singular value < tol
% are considered to lie in the null space.
% If no tolerance value is given, the default
%   tol = max(m,n) * max(s) * eps(class(A));
% is assumed.

[m,n] = size(A);

[~,S,V] = svd(A,0);
if m > 1, s = diag(S);
elseif m == 1, s = S(1);
else s = 0;
end

if nargin < 2
  tol = max(m,n) * max(s) * eps(class(A));
end
r = sum(s > tol);
Z = V(:,r+1:n);
end