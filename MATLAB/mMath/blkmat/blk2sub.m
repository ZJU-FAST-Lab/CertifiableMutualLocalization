function sub = blk2sub(blocks, block_sizes)
% BLK2SUB Return a vector of subscripts corresponding to the specified blocks.
% sub = blk2sub(blocks, block_sizes)
%
% e.g., blk2sub([2 5], [2 1 2 1 2]) = [3 7 8].

blocks = blocks(:)';

if numel(block_sizes)==1
  % If the block size is regular, expand block_sizes as much as necessary
  % according to input block indeces
  block_sizes = repmat(block_sizes,1,max(blocks));
end
block_sizes = block_sizes(:)';
skip = [0 cumsum(block_sizes)];
start = skip(blocks)+1;
fin = start + block_sizes(blocks) - 1;
sub = [];
for j=1:length(blocks)
  sub = [sub start(j):fin(j)];
end
end