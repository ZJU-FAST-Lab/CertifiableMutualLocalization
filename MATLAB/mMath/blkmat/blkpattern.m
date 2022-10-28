classdef blkpattern
  %blkpattern Class that store the metadata in a block-matrix
  %   This class handles all the metadata in relation to a block-matrix.
  
  properties
    % Metadata
    sizes
    dict
    % Flags, set once at the constructor and never modified again
    is_regular
    is_labeled
  end
  properties (Dependent)
    labels
  end
  
  methods
    function this = blkpattern(varargin)
      if nargin == 0
        this.sizes = [];
        this.dict = struct();
        
      elseif nargin == 1
        if isa(varargin{1}, 'blkpattern')
          % Copy constructor
          this = varargin{1};
          return % When clone copy is done, we can skip subsequent steps
        end
        
      else
        % There are more than 1 input argument
        dims = varargin{1}; sizes = varargin{2};
        [this.sizes,this.dict] = setupStructure( dims,sizes );
      end
      
      % Set regularity flag
      this.is_regular = (numel(unique(this.sizes))==1);
      % Set flag for labelled pattern
      this.is_labeled = ~isempty(fieldnames(this.dict));
    end
       
    % Arithmetic
    function is = eq(a,b)
      % Compare the two different blk patterns. If the patterns are
      % labeled, the labeling must be consistent for equality.
      is = all(a.sizes == b.sizes) && all(a.labels == b.labels);
    end
    
    % Utilities
    function l = get.labels(this)
      l = cell2mat(fieldnames(this.dict));
    end
    
    function intIdxs = numSubs( this, idxs )
      if isnumeric( idxs )
        intIdxs = idxs;
      else
        if strcmp(idxs,':')
          intIdxs = 1:numel(this.sizes);
        else
          % Check all the labels are contained in the dictionary
          labelsNotInDict = setdiff(idxs,this.labels);
          assert(isempty(labelsNotInDict),...
            'Labels %s don''t exist in the matrix dimension',labelsNotInDict);
          % Read indices corresponding to labels from the dictionary
          n = numel(idxs);
          intIdxs = zeros(1,n);
          for j=1:n
            l = idxs(j);
            intIdxs(j) = this.dict.(l);
          end
        end
      end
    end
    
    function idxs = block( this, blkIdxs )
      idxs = blk2sub(blkIdxs, this.sizes);
    end
    
  end
  
end

function [sizes,dict] = setupStructure( arg1, sizes )
% Preallocate dictionary
dict = struct();

if ischar(arg1)
  rtags = arg1;
  nblks = numel(rtags);
  for i = 1:nblks
    l = rtags(i);
    dict.(l) = i;
  end
else
  nblks = arg1;
  if isempty(nblks)
    % If empty nblks, this is take from dimension of rsizes
    nblks = numel(sizes);
  end
end
if isscalar(sizes)
  % If we are given a single block-size, repeat this
  sizes = repmat(sizes, 1, nblks);
end
assert(nblks==numel(sizes),...
  'Number of blocks (%d) and list of block dims (%d) should be consistent',...
  nblks, numel(sizes));
end