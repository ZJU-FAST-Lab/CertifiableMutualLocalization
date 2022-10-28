function Q = orth(A,tol)
% Q = ORTH(A)
% Q = ORTH(A,tol)
% 
% Q = ORTH(A) is an orthonormal basis for the range of A obtained
% from the singular value decomposition.
% Namely, all those right singular vectors with singular value > tol
% are considered to lie in the null space.
% If no tolerance value is given, the default
%   tol = max(size(A)) * S(1) * eps(class(A));
% is assumed.

[Q,S] = svd(A,'econ'); %S is always square.
if ~isempty(S)
    S = diag(S);
    if nargin < 2
      tol = max(size(A)) * S(1) * eps(class(A));
    end
    r = sum(S > tol);
    Q = Q(:,1:r);
end