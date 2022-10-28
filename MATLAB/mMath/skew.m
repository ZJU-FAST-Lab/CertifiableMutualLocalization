function M = skew(v)
% M = SKEW(v)
% Creates a skew symmetrix matrix from a 3-vector as
%   M = [ 0  -v3 +v2
%        +v3  0  -v1
%        -v2 +v1  0 ]
% Some overloaded cases are:
% - v is a cell (list) of three elements (with consistent dimensions)
% - v is a 3-block-vector (regular)
% 
%   See also INV_SKEW, BLKMAT.

% This only makes sense for 3-dimensional vectors or lists
assert(numel(v)==3)

switch class(v)
  case 'cell'
    % The input is a 3-element cell list
    s  = size(v{1});
    Os = zeros(s);
    M  = [ Os   -v{3}  v{2};
           v{3}  Os   -v{1};
          -v{2}  v{1}   Os];
  case 'blkmat'
    % The input is a 3-element block-vector
    assert(isregular(v),'Blk-matrix must be regular');
    
    Os = zeros(blksize(v));
    M = [  Os   -v(3)  +v(2)
         +v(3)   Os    -v(1)
         -v(2)  +v(1)   Os  ];
    
    M = blkmat(3,3,rowsize(v),colsize(v),M);
  otherwise
    % The usual naive case: v\in\Re^3
    M = [  0    -v(3)  +v(2)
         +v(3)    0    -v(1)
         -v(2)  +v(1)    0  ];
end

end
