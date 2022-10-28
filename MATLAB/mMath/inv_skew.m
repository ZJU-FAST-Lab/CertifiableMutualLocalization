function v = inv_skew(M,nocheck)
% v = INV_SKEW(M)
% v = INV_SKEW(M,'nocheck')
% For a skew-symmetric matrix extracts v so that skew(v)=M
% Skew-symmetry is checked internally unless 'nocheck' option is given
% 
% See also SKEW, BLKMAT.

if nargin == 1
  % Check that the input matrix is skew-symmetric
  tol = 1e-6;
  assert(norm(M(:,:)+M(:,:)','fro')<tol,'Matrix must be skew-symmetric');
else
  % Do not check skew-symmetry (assume skew-symmetry is assured externally)
  if ~strcmp(nocheck,'nocheck')
    error('Invalid 2nd argument, did you mean ''nocheck''?');
  end
end

switch class(M)
  case 'blkmat'
    % The input is a 3-element block-vector
    assert(isregular(M),'Blk-matrix must be regular');   
    v = blkmat(3,1,rowsize(M),colsize(M));
  otherwise
    v = zeros(3,1);
end
% Extract values
v(1) = -M(2,3);
v(2) = +M(1,3);
v(3) = -M(1,2);

end