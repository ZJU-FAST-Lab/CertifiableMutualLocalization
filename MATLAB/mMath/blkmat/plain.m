function M = plain(A)
% M = plain(M)
% Returns a usual Matlab-matrix object
% This can be used to have a uniform interface in blkmat operators.
% 
% Note that we overload this method for blkmat-like classes.
% Thus, if we call plain(M) where M is an usual Matlab variable,
% this method is called. If instead we call plain(M) with M an object
% of a specific blkmat-like class, the custom plain() method is entered.
% 
% See also blkmat, graphMat.

if ismatrix(A)
  M = A;
else
  error('The input for plain must be a matrix-like object');
end
end