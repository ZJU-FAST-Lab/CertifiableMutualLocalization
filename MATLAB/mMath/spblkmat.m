function M = spblkmat(nrows,ncols,rsize,csize,nzmax)
% M = SPBLKMAT(nrows,ncols,rsize,csize,nzmax)
% SPBLKMAT Create BLKMAT with sparse internal storage.
% Allocate space to eventually hold NZMAX non-empty (full) blocks.
% 
% See also BLKMAT, SPALLOC.

numelBlock = rsize*csize;
internalM = spalloc( nrows*rsize, ncols*csize, numelBlock*nzmax );
M = blkmat(nrows,ncols,rsize,csize,internalM);
  
end