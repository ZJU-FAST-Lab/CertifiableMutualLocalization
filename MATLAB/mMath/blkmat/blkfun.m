function A = blkfun(fun, B)
% A = blkfun(fun, B)
% Apply a function to each block of a block matrix.
%
% See also: cellfun, ARRAYFUN, STRUCTFUN, SPFUN

assert(isa(B,'blkmat'));
A = zeros(nrows(B),ncols(B));
for i=1:nrows(B)
  for j=1:ncols(B)
    A(i,j) = fun(B(i,j));
  end
end