function R = rot(n)
% Generate random rotation in SO(n), that is, R'R=I and det(R)=+1
%   R = rot()
% By default n=3

if nargin==0
  n=3;
end
[R,~] = qr(randn(n));
% Assert determinant is +1
if det(R) < 0
  R(:, [1 2]) = R(:, [2 1]);
end