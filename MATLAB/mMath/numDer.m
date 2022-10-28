function D = numDer( F, X0, h, N )
% D = numDer( F, X0 )
% D = numDer( F, X0, h )
% 
% Compute numerical (tensor) differentiation
% By tensor differentiation we mean keeping the same number of indeces
% to refer to the output elements as the sum of indeces in the inputs:
% If we have a function
%   F : R^{a1,...,an} -> R^{b1,...,bm}
% whose input is a n-dim array and whose output is a m-dim array,
% the derivative is a (n+m)-dim array of dimensions:
%   D is R^{b1,...,bm,a1,...,an}
% 
% The singleton dimensions (dimensions of size 1) are ignored
% to produce more compact arrays.
% 
% Special cases -> From its definition, this operator falls back to:
%   - Scalar derivative if both F and X0 are scalar
%   - Gradient (as column vector) if F is scalar and X0 is vector
%   - Jacobian (f's are rows and x's are columns) if F and X0 are vectors
% 
% 
% D = numDer( F, X0, h, N )
% 
% Compute N-th order numerical (tensor) differentiation
% That means, if we have a function
%   F : R^{a1,...,an} -> R^{b1,...,bm}
% whose input is a n-dim array A = {a1,...,an}
% and whose output is a m-dim array B = {b1,...,bm},
% the derivative is a (m+n^N)-dim array of dimensions:
%   D is R^{B,A^N} = R^{B,A,...,A}
% 
% Special cases -> If N=2, this operator falls back to:
%   - Hessian if F is scalar and X0 is vector

if nargin < 3
  h = 1e-5;
end
if nargin < 4
  N = 1;
end

% Higher-order derivatives can be defined in terms of lower derivatives
if N > 1
  % Define auxiliary inner function (a derivative of order N-1)
  derF = @(X)numDer( F, X, h, N-1 );
  % Recursive call on the auxiliary function
  D = numDer( derF, X0, h );
  return
end

% Get dimensions of input and output
F0 = F(X0);
sizeF = size(F0);
sizeX = size(X0);

% Store intermediate result as vectorized
Dvec = zeros( numel(F0), numel(X0) );
% Compute scalar-wise derivatives for each input variable
for i = 1:numel(X0)
  % Define unitary vector in the i-th (vectorized) dimension
  e = zeros(sizeX); e(i)=1;
  % Apply central difference rule
  DF_Xi = (1/h)*0.5*(F(X0+h*e) - F(X0-h*e));
  % Save vectorized version of the derivative
  Dvec(:,i) = DF_Xi(:);
end

% Reshape computed results to original dimensions
Dsize = [sizeF(sizeF~=1),sizeX(sizeX~=1)];
if numel(Dsize)<2 % Derivative is a scalar or a vector
  D = Dvec(:);
else
  D = reshape(Dvec,Dsize);
end
end