function X = avec(x,varargin)

% TODO: add dimension argument
if nargin == 1
  % assume square 2D matrix
  dim = sqrt(numel(x));
  assert( dim-fix(dim) == 0, 'Is matrix square?' );
  dim = [dim,dim];
else
  dim = cell2mat(varargin);
end
assert( all( dim-fix(dim) == 0 ) )
X = reshape(x,dim);

end