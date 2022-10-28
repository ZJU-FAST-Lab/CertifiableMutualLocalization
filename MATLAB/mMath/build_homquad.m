function Q_ = build_homquad( Q,b,c )
% Q_ = build_homquad( Q,b,c )
% Build homogenous quadratic form from quadratic, linear and independent
% terms in a non-quadratic form.
% The symmetric representation is taken by symmetrization
% 
% See also makehom,makenonhom.

assert(issquare(Q))
n = size(Q,1);

Q_ = [Q,b;zeros(1,n),c];
Q_ = symmetrize(Q_);