function Pmn = build_vecProj( m,n )
% Pmn = build_vecProj( m,n )
% 
% Builds the projection matrix that fulfills the relation
%   Pmn*vec(A) = vec(A')
% where A is a mxn matrix

if nargin == 1
  n = m;
end

permIdxs = vec(reshape(1:(m*n),m,n)');

Pmn = eye(m*n);
Pmn = Pmn(:,permIdxs);