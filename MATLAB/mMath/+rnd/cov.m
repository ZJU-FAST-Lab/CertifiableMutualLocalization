function A = cov(n,nze)
% Generate random nxn symmetric positive semidefinite matrix
%   A = cov(n)
% This is useful for getting covariance matrices

[U,S,~] = svd(randn(n));
if nargin > 1
  assert(nze <= n,'Number of zero eig must be <= number of eig')
  s = diag(S);
  s = [ s(1:n-nze); zeros(nze,1) ];
  S = diag(s);
end
A = U*S*U';

% To avoid numerical issues, the result is symmetrized!
% Otherwise epsilon magnitude errors appear in norm(A-A')
A = symmetrize(A);