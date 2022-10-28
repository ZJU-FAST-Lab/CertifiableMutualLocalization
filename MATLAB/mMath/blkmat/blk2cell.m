function c = blk2cell(B)
%BLK2CELL Convert block numeric array into cell array.
% TODO: Extend options to the same as num2cell
%
% See also num2cell.

% preallocate cell array
c = cell(nrows(B),ncols(B));

for i=1:nrows(B)
  for j=1:ncols(B)
    c{i,j} = B(i,j);
  end
end

end