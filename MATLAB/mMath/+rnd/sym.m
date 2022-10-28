function X = sym(n)
% Generate random nxn symmetric matrix
%   X = sym(n)

X = symmetrize( randn(n) );