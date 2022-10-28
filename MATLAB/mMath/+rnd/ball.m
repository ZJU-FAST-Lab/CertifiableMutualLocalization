function x = ball(d,n,R)
% Points in the interior of a sphere in R^d
% x = ball(d,n,R)
% Parameters:
%   d - dimension of the ball: x \in Re^d
%   n - number of points to generate
%   R - radius of the ball: ||x||_2 <= R
% 
% This function can be useful to generate bounded spherical noise.
% 
% Note we take a random direction and then sample a point uniformly
% between the center and the boundary of the ball in that direction.
% As a result, points will tend to be more concentrated in the center.

% directions (from normal distribution)
dir  = snormalize(randn(d,n));
% distances between 0 and R (uniform)
dist = R*rand(1,n);

% compose points
x = repmat(dist,d,1).*dir;

end