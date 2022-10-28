function X = vec_inv(x)

% TODO: add dimension argument
if nargin == 1
  dim = sqrt(numel(x));
  assert( dim-fix(dim) == 0 ); % Integer value
  X = reshape(x,dim,dim);
end

end