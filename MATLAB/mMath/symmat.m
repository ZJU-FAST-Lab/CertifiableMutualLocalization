function x = symmat(var,dims,varargin)
% x = SYMMAT(var,n)
% Create a symbolic multidimensional array of basename var
% with dims(1)x...xdims(end) elements
% 
% SYMMAT(___, set) sets the assumption that the variables belong to a set.
% Here, set can be 'real', 'positive', 'integer', or 'rational'.
% 
% See also SYM, SYMVEC

numDims = numel(dims);
if any(dims>=10)
  varsub = repmat('_%d',1,numDims);
else
  varsub = ['_',repmat('%d',1,numDims)];
end
varname = [var,varsub];
x = sym(varname,dims,varargin{:});