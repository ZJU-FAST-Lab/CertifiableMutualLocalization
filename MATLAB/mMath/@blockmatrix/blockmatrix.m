function M = blockmatrix(nrows, ncols, rsize, csize)
% BLOCKMATRIX Constructor for the blockmatrix class.
%
% M = BLOCKMATRIX(NROWS, NCOLS, RSIZE, CSIZE) creates an NROWS by NCOLS blockmatrix,
% where the size of each block is RSIZE by CSIZE. This is called a regular blockmatrix.
%
% It is also possible for the blocks to have different sizes along a dimension.
% M = BLOCKMATRIX([], NCOLS, RSIZES, CSIZE) creates a NROWS by NCOLS blockmatrix,
% where NROWS = LENGTH(RSIZES), and the size of block (i,j) is RSIZES(I) by CSIZE.
% This is called a row irregular, column regular blockmatrix.
% M = BLOCKMATRIX(NROWS, [], RSIZE, CSIZES) creates a NROWS by NCOLS blockmatrix,
% where NCOLS = LENGTH(CSIZES), and the size of block (i,j) is RSIZE by CSIZES(J)
% This is called a row regular, column irregular blockmatrix.
% M = BLOCKMATRIX([], [], RSIZES, CSIZES) creates a NROWS by NCOLS blockmatrix,
% where the size of block (i,j) is RSIZES(I) by CSIZES(J).
% This is called a fully irregular matrix.
%
% Note that M = BLOCKMATRIX(NROWS, NCOLS, RSIZE, CSIZE) is equivalent to
% M = BLOCKMATRIX([], [], RSIZE*ONES(1,NROWS), CSIZE*ONES(1,NCOLS)),
% except that the latter is considered irregular.
%
% Use 'x = M(i,j)' to assign block (i,j) to x.
% If this block does not exist, an error will be raised.
% If ncols(M)=1, it suffices to write 'x = M(i)', and similarly for row vectors.
% You can also use ranges, e.g., M(:,:) or M(:) or M(1:3)
%
% Use 'M(i,j) = x' to assign x to block (i,j).
% If this block does not exist, it will be created,
% providing M is regular in the out-of-bounds dimension.
% (If M is irregular, an error will be raised.)
%
% The following functions are defined in the obvious way:
% NROWS(M), NCOLS(M)
% ROWSIZES(M), COLSIZES(M) if irregular
% ROWSIZE(M), COLSIZE(M) if regular
%
% No other operations can be performed on blockmatrices.
% To do arithmetic, you need to operate on the underlying scalar matrix,
% which can be accessed using M(:,:).
%
% EXAMPLES
%   M = blockmatrix(2, 2, 2, 2)
%   M(1,1) = eye(2); M(2,2) = -eye(2).
% 
%     1     0     0     0
%     0     1     0     0
%     0     0    -1     0
%     0     0     0    -1
%
%   M = blockmatrix([], [], [2 3], [2 1]) 
%   M(:,:) = rand(5, 3) looks like this:
%
%    x x   x
%    x x   x
%
%    x x   x
%    x x   x
%    x x   x

% Written by Kevin Murphy (www.cs.berkeley.edu/~murphyk), 21 October 1999


if nargin == 0 % default constructor
  M.rsizes = [];
  M.csizes = [];
  M.row_regular = 0;
  M.col_regular = 0;
  M.M = [];
  M = class(M, 'blockmatrix');
  
elseif isa(nrows, 'blockmatrix')
  M = nrows; % identity function
  
else
  if isempty(nrows)
    M.rsizes = rsize;
    M.row_regular = 0;
  else
    M.rsizes = repmat(rsize, 1, nrows);
    M.row_regular = 1;
  end
  
  if isempty(ncols)
    M.csizes = csize;
    M.col_regular = 0;
  else
    M.csizes = repmat(csize, 1, ncols);
    M.col_regular = 1;
  end
  M.M = zeros(sum(M.rsizes), sum(M.csizes));
  M = class(M, 'blockmatrix');
end
