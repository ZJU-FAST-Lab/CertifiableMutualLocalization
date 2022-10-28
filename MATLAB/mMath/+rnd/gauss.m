function x = gauss(mu,cov,numSamples)
% x = gauss(mu,cov)
% Generate values from a bivariate normal distribution with
% specified mean vector and covariance matrix.
% See also Example 2 in randn.
% 
% X = gauss(mu,cov,numSamples)
% Same as above but produces numSamples as columns of X
% 
% See also randn.

assert(isvector(mu),'The mean must be a vector')
assert(issymmetric(cov),'Covariance must be symmetric');
U = chol(cov);
n = numel(mu);
if nargin < 3
  x = mu + U'*randn(n,1);
else
  x = repmat(mu,1,numSamples) + U'*randn(n,numSamples);
end