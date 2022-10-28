% blkmat Blockmatrix class.
% This class works as a convenient wrapper to improve code readability
% when working with matrices that are conceptually distributed blockwise.
% The main purpose of the created object is to store a usual matrix while
% keeping a blockwise interface, where subscript indices refer to blocks.
%
% The indexing (subsasgn,subsref) methods are modified so that the objects
% behaves externally as a matrix:
% 
%   - Read block
%   Use `x = M(i,j)` to assign block (i,j) to x.
%   If this block does not exist, an error will be raised.
%   If ncols(M)=1, it suffices to write 'x = M(i)', and similarly for row vectors.
%   You can also use ranges, e.g., M(:,:) or M(:) or M(1:3) or M(2:end)
%
%   - Write block
%   Use `M(i,j) = x` to assign x to block (i,j).
%   If this block does not exist, an error will be raised.
% 
% We enrich the object with the ability to tag the blocks with letters, so
% that letter indices and struct syntax can be used to access the blocks
% (as long as the blkmat is built as labeled, see constructors):
% 
%   - Access block with letter indices
%   Use `M('a','mn')` to access the block tagged *a* in the row dimension
%   and *b* in the second dimension. 
%   The struct syntax `M.am` is equivalent to `M('a','m')`. This syntax
%   allows only for single-index subscripts (not lists of indices).
% 
% 
% CONSTRUCTORS
% ------------
% 
%   M = blkmat(NROWS, NCOLS, RSIZE, CSIZE)
% Create an NROWS by NCOLS blockmatrix, where the size of each block
% is RSIZE by CSIZE. This is called a regular blockmatrix.
%
% It is also possible for the blocks to have different sizes along a dimension.
% In this case we give a vector of block-sizes:
%   M = blkmat([], NCOLS, RSIZES, CSIZE)
%   M = blkmat(NROWS, [], RSIZE, CSIZES)
%   M = blkmat([], [], RSIZES, CSIZES)
% If the number of blocks is given in addition to the vector of sizes,
% we keep the same behavior (with additional checking of the consistency).
% If the block-sizes in a certain dimension are not all the same,
% this dimension is called non-regular.
% 
%   M = blkmat(RLABELS, CLABELS, ...)
% Create a *labeled* blkmat, where the blocks can be referred to by its
% tag in addition to their numeric index. See examples.
% 
% If the matrix structure is symmetric, we can give the dimensions and/or
% sizes only once:
%   M = blkmat(NBLOCKS, SIZE)
%   M = blkmat([], SIZES)
%   M = blkmat(LABELS, SIZES)
% 
% The internal value of the blkmat can be initialized giving a plain matrix
% (with consistent dimensions) as the last argument:
%   M = blkmat( ..., M0 )
%
% Other useful constructors are:
%   B = blkmat(A)
%   Clone the object A (same content).
% 
%   B = blkmat(A,M0)
%   Copy the structure of A but stores data from M0.
% 
%   B = blkmat(RPATTERN,CPATTERN)
%   B = blkmat(RPATTERN,CPATTERN,M0)
%   Create a matrix with the input patterns (objects of class blkpattern).
% 
% USEFUL METHODS
% --------------
% 
% The following functions are defined in the obvious way:
% nrows(M), ncols(M)
% rowsizes(M), colsizes(M)
% rowsize(M), colsize(M) (only useful if regular)
%
% 
% OPERATORS
% ---------
% Usual matrix operations are overloaded for blkmat checking for internal
% consistency of the block sizes and, in general, returning the
% corresponding new blkmat object.
%
% EXAMPLES
%   M = blkmat(2, 2, 2, 2)
%   M(1,1) = eye(2); M(2,2) = -eye(2).
%
%     1     0     0     0
%     0     1     0     0
%     0     0    -1     0
%     0     0     0    -1
%
%   M = blkmat([], [], [2 3], [2 1])
%   M(:,:) = rand(5, 3) looks like this:
%
%    x x   x
%    x x   x
%
%    x x   x
%    x x   x
%    x x   x
% 
%   M = blkmat('AB','CD', 3, 3)
%   M.AC; M('A','C'); M('AB','C'); M(:,'D')
% 
%        C       D
%      -----   -----
%        
%    | x x x   x x x
%   A| x x x   x x x
%    | x x x   x x x
%
%    | x x x   x x x
%   B| x x x   x x x
%    | x x x   x x x
%
% Notes:
% - The usual functionalities for matrix objects are extended.
%   For avoiding confusion, some of the traditional methods are
%   prepended the *blk* word, e.g. blknumel().
% 
% See also: blkpattern.

% Inspired by Kevin Murphy (www.cs.berkeley.edu/~murphyk), 21 October 1999
% This version by Jesus Briales, 3 May 2016

classdef blkmat
  
  properties
    rpattern, cpattern
    storage
  end
  
  methods   
    function this = blkmat( varargin )
      % (nrows, ncols, rsize, csize, M)
      % blkmat
      % Possible inputs are:
      % - NROWS, NCOLS, RSIZE, CSIZE
      % - [], [], RSIZES, CSIZES
      
      if nargin == 0
        % No arguments, default constructor
        this.rpattern = blkpattern();
        this.cpattern = blkpattern();
        this.storage = [];
        
      elseif nargin == 1
        % One single input argument
        if isa(varargin{1}, 'blkmat')
          % Copy constructor
          this = varargin{1};
          return % When clone copy is done, we can skip subsequent steps
        else
          % Build single block matrix from numeric input
          M = varargin{1};
          assert(isnumeric(M),'The input should be a numeric matrix');
          s = size(M);
          this.rpattern = blkpattern(1,s(1));
          this.cpattern = blkpattern(1,s(2));
          this.storage = M;
        end
        
      else
        % There are more than 1 input argument
        if isa(varargin{1}, 'blkmat')
          % If 1st argument is blkmat, copy the input structure
          this = varargin{1};
          % There is extra argument giving content initialization
          M = varargin{2};
          
        elseif isa(varargin{1}, 'blkpattern')
          % Directly give blkpatterns as inputs
          this.rpattern = varargin{1};
          if nargin >= 2 && isa(varargin{2}, 'blkpattern')
            this.cpattern = varargin{2};
          else
            this.cpattern = varargin{1}; % Symmetric pattern
          end
          if ~isa(varargin{end}, 'blkpattern')
            % If the last argument is NOT blkpattern, it must be the init.
            M = varargin{end};
          else
            M = 0;
          end
          
        else
          if nargin <= 3
            % The symmetric case, same dims and blk-sizes for rows and columns
            dims = varargin{1}; bsizes = varargin{2};
            pattern = blkpattern(dims,bsizes);
            this.rpattern = pattern; this.cpattern = pattern;
          else
            % Set row structure
            rdims = varargin{1}; rsizes = varargin{3};
            this.rpattern = blkpattern(rdims,rsizes);
            % Set col structure
            cdims = varargin{2}; csizes = varargin{4};
            this.cpattern = blkpattern(cdims,csizes);
          end
          
          % The extra argument for initialization can be the 3rd or 5th,
          % so it is the last if there is an odd number of arguments
          if mod(nargin,2)
            % If odd, there is extra argument giving content initialization
            M = varargin{end};
            % If the input is already blkmat
            if isa(M,'blkmat')
              % Assert the blk dimensions are the same
              assert( all(rowsizes(this)==rowsizes(M)) && ...
                      all(colsizes(this)==colsizes(M)) );
              M = plain(M);
            end
          else
            % If even, no initialization, set to 0
            M = 0;
          end
        end
        
        % Initialize storage field
        assert(isscalar(M) || all(size(M)==size(this)),...
          'blkmat: Wrong dim of the initialization matrix')
        if isscalar(M)
          this.storage = M*ones(size(this));
        else
          this.storage = M;
        end
      end
    end
       
    % Set the number of arguments in subscript operations
    % Avoid conflicts with (overloaded) numel method
    function n = numArgumentsFromSubscript(this,~,callingContext)
      switch callingContext
        case matlab.mixin.util.IndexingContext.Statement
          n = 1; % nargout for indexed reference used as statement
        case matlab.mixin.util.IndexingContext.Expression
          n = 1; % nargout for indexed reference used as function argument
        case matlab.mixin.util.IndexingContext.Assignment
          n = 1; % nargin for indexed assignment
      end
    end
    
    function ind = end(this,k,n)
      if k==1
        ind = nrows(this);
      elseif k==2
        ind = ncols(this);
      else
        error('Only 2 subscripts are used by blkmat');
      end
    end
    
    function [rows,cols] = extract_indices( this, S )
      if S.type == '.'
        % If field accessor, convert into usual indexing
        % (only one-char index per dimension can be used)
        S.type = '()';
        S.subs = num2cell(S.subs);
      end
      if strcmp(S.type,'()')
        % Complete indexing for matrix
        switch numel(S.subs)
          case 1
            % Case a single index is given (blk-vector)
            if ncols(this)==1 % col vector
              S.subs = {S.subs{1}, 1};
            elseif nrows(this) == 1 % row vector
              S.subs = {1, S.subs{1}};
            else
              error('This is not a blk-vector, specify two indices for matrices');
            end
          case 2
            % Everything fine
          otherwise
            error('Blk-mat has maximum 2 dimensions');
        end
        
        % Convert all idxs to numeric idxs
        % For matrix, there are exactly 2 indices (row and column)
        rows = this.rpattern.numSubs( S.subs{1} );
        cols = this.cpattern.numSubs( S.subs{2} );
        
      end
    end
    
    function this = subsasgn(this,S,B)
      
      if length(S)==1 % A(i,j), A.ab, A('a','b')
        assert(strcmp(S.type,'()') || strcmp(S.type,'.'),...
          'Unrecognized subscript operator %s', S.type);
        % Get numeric indices from subscript
        [rows, cols] = extract_indices(this, S);
        % Check out-of-bound indeces
        if max(rows) > nrows(this) || max(cols) > ncols(this)
          error(['index ' num2str(rows) ', ' num2str(cols) ' is out of bounds']);
        end
        % Set internal storage from input matrix
        this.storage(this.rpattern.block(rows),...
                     this.cpattern.block(cols)) = B;
      else
        error(['Unrecognized subscript operator ' S]);
      end
    end
    
    function varargout = subsref(this,S)
      
      if length(S)==1 % A(i,j), A.ab, A('a','b')
        assert(strcmp(S.type,'()') || strcmp(S.type,'.'),...
          'Unrecognized subscript operator %s', S.type);
        % Get numeric indices from subscript
        [rows, cols] = extract_indices(this, S);
        % Check out-of-bound indeces
        if max(rows) > nrows(this) || max(cols) > ncols(this)
          error(['index ' num2str(rows) ', ' num2str(cols) ' is out of bounds']);
        end
        % Set output matrix from blocks
        varargout{1} = this.storage(this.rpattern.block(rows),...
                                    this.cpattern.block(cols));
      else
        error(['Unrecognized subscript operator ' S]);
      end
    end
        
    % Unary operators
    function B = uminus(A)
      B = A;
      B.storage = -A.storage;
    end
    function At = ctranspose(A)
      At = blkmat([],[],colsizes(A),rowsizes(A),A.storage');
    end
    function t = trace(A)
      t = trace(A.storage);
    end
    
    % Binary operators
    function [rpat,cpat] = sumPattern(A,B)
      % Valid cases are:
      % - A and B blkmat, with compatible size and blksize
      assert( A.rpattern == B.rpattern && A.cpattern == B.cpattern )
      rpat = A.rpattern; cpat = A.cpattern;
    end
    
    function [rpat,cpat] = prodPattern(A,B)
      % Valid cases are:
      % - A and B blkmat, with compatible size and blksize
      % - A or B scalar
      if isscalar(plain(A))
        rpat = B.rpattern; cpat = B.cpattern;
      elseif isscalar(plain(B))
        rpat = A.rpattern; cpat = A.cpattern;
      else
        assert( A.cpattern == B.rpattern )
        rpat = A.rpattern; cpat = B.cpattern;
      end
    end
    
    function C = plus(A,B)
      
      % Force blkmat matrices (in case one of the inputs is numeric)
      A = blkmat(A); B = blkmat(B);
      
      % Compute result structure
      [rpat,cpat] = sumPattern(A,B);
      
      % Sum of internal storages
      temp = plain(A) + plain(B);
      
      % Store result in blkmat with the same structure
      C = blkmat(rpat,cpat,temp);
    end
    
    function C = minus(A,B)
      
      % Force blkmat matrices (in case one of the inputs is numeric)
      A = blkmat(A); B = blkmat(B);
      
      % Compute result structure
      [rpat,cpat] = sumPattern(A,B);
      
      % Difference of internal storages
      temp = plain(A) - plain(B);
      
      % Store result in blkmat with the same structure
      C = blkmat(rpat,cpat,temp);
    end  
        
    function C = mtimes(A,B)
      % C = mtimes(A,B)
      % Some semantics in the product with block matrix:
      % - Product with a scalar simply scales all the content
      % - If the result is a single block (scalar or matrix)
      %   return raw Matlab matrix
      
      % Force blkmat matrices (in case one of the inputs is numeric)
      A = blkmat(A); B = blkmat(B);
           
      % Compute result structure
      [rpat,cpat] = prodPattern(A,B);
      
      % Product of internal storages
      temp = plain(A)*plain(B);
      
      % Store result in blkmat with the resulting structure
      C = blkmat(rpat,cpat,temp);
      if numel(C) == 1
        % Result is a single block, return as simpler numeric matrix
        C = plain(C);
      end
    end
    
    function C = mldivide(A,B)
      % C = mldivide(A,B)
      
      % Force blkmat matrices (in case one of the inputs is numeric)
      A = blkmat(A); B = blkmat(B);
      
      % Compute result structure
      [rpat,cpat] = prodPattern(A,B);
      
      % Left division of internal storages
      temp = plain(A)\plain(B);
      
      % Store result in blkmat with the resulting structure
      C = blkmat(rpat,cpat,temp);
      if numel(C) == 1
        % Result is a single block, return as simpler numeric matrix
        C = plain(C);
      end
    end
    
    function C = mrdivide(A,B)
      % C = mldivide(A,B)
      
      % Force blkmat matrices (in case one of the inputs is numeric)
      A = blkmat(A); B = blkmat(B);
            
      % Compute result structure
      [rpat,cpat] = prodPattern(A,B);
      
      % Right division of internal storages
      temp = plain(A)/plain(B);
      
      % Store result in blkmat with the resulting structure
      C = blkmat(rpat,cpat,temp);
      if numel(C) == 1
        % Result is a single block, return as simpler numeric matrix
        C = plain(C);
      end
    end
    
    function r = rowsize(A),  r = A.rpattern.sizes(1); end
    function r = rowsizes(A), r = A.rpattern.sizes; end
    function r = nrows(A),    r = length(A.rpattern.sizes); end
    function c = colsize(A),  c = A.cpattern.sizes(1); end
    function c = colsizes(A), c = A.cpattern.sizes; end
    function c = ncols(A),    c = length(A.cpattern.sizes); end
    function s = size(A)
      % TODO: Refactor so that size counts nrows and ncols
      % For complete size use size(A(:,:)) instead
      s = [sum(A.rpattern.sizes), sum(A.cpattern.sizes)];
    end
    function n = numel(A), n = nrows(A)*ncols(A); end
    function s = blksize(A)
      % NOTE: regular blkmat should be asserted externally
      s = [rowsize(A),colsize(A)];
    end
    function is = isregular(A), is = A.rpattern.is_regular && A.cpattern.is_regular; end
    function l = labels(this)
      l = {cell2mat(fieldnames(this.rpattern.dict)),...
           cell2mat(fieldnames(this.cpattern.dict))};
    end
    
    function M = plain(this), M = this.storage; end

    function is = isempty(A)
      is = isempty(A.storage);
    end
    function E = blkisempty(A)
      % E = blkisempty(A)
      % Check for each block if all elements are zero
      nr = nrows(A); nc = ncols(A);
      E = false(nr,nc);
      for i=1:nr
        for j=1:nc
          B = A.storage(blk2sub(i,rowsizes(A)), blk2sub(j,colsizes(A)));
          E(i,j) = all(all(B==0));
        end
      end
    end
    
    function disp(this)
      % Display the metadata about blk-matrix structure only
      % To display the matrix content, use plain method
      
      c = class(this.storage);
      if this.rpattern.is_regular && this.cpattern.is_regular
        s = sprintf('%dx%d %s-matrix of %dx%d blocks',...
          nrows(this), ncols(this), c, rowsize(this), colsize(this));
      elseif this.rpattern.is_regular && ~this.cpattern.is_regular
        s = sprintf('%dx- %s-matrix of %d x [%s] blocks',...
          nrows(this), c, rowsize(this), num2str(colsizes(this)));
      elseif ~this.rpattern.is_regular && this.cpattern.is_regular
        s = sprintf('-x%d %s-matrix of [%s] x %d blocks',...
          ncols(this), c, num2str(rowsizes(this)), colsize(this));
      else
        s = sprintf('-x- %s-matrix of [%s] x [%s] blocks',...
          c, num2str(rowsizes(this)), num2str(colsizes(this)));
      end
      
      % Add information about labels, if the blkmat is labeled
      strLabels = [];
      if this.cpattern.is_labeled || this.rpattern.is_labeled
        strLabels = sprintf(', labeled as [%s] x [%s]',...
          num2str(cell2mat(fieldnames(this.rpattern.dict))),...
          num2str(cell2mat(fieldnames(this.cpattern.dict))) );
      end
      
      % Print to screen
      fprintf('\t%s%s\n',s,strLabels);
    end
    
    function spy(A)
      spy(A.storage);
      % If exists spyblck, add lines between blocks
      % TODO: Check function exists
      % TODO: Should implement spyblk in the blkpattern class
      plt.spyblk(rowsizes(A),colsizes(A));
      
      % Adjust borders
      a = axis;
      axis(a + [+0.5 -0.5 +0.5 -0.5]);
    end
    
    function out = full(A)
      out = full(A.storage);
    end
    
  end

end