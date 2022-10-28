function main( )

% Tests for blkmat functionalities

%% Constructors
% Regular matrix with
% - number of rows and columns
% - (homogeneous) size of blocks in each dimension
B = blkmat(2,3,4,5);
B = blkmat(2,3,4,5,randn(2*4,3*5));

% Non-regular matrix
B = blkmat([],[],[2 3], [4 5 6]);
% Non-regular matrix with redundant dims
B = blkmat(2,3,[2 3], [4 5 6]);

% Non-regular matrix with symmetric tagged dimensions
B = blkmat('QT',[2 3]);
B = blkmat('QT',[2 3],randn(5,5));

% Single block matrix from numeric matrix
B = blkmat(randn(3,4));

% Clone matrix
B = blkmat(2,3,4,5);
C = blkmat(B);

% Clone structure with different content
B = blkmat(2,3,4,5);
C = blkmat(B,randn(2*4,3*5));

% Regular matrix with non-symmetric labels
B = blkmat('ab','cd',2,3);
B = blkmat('ab','cd',2,3,randn(2*2,2*3));

%% Operators
% Sum of block-matrices
A = blkmat(2,2,2,2,randn(4,4));
B = blkmat(2,2,2,2,randn(4,4));
A + B;
% Sum of 

% Product of non-regular matrices
A = blkmat([],[],[1 2],[3 4 5]);
B = blkmat([],[],[3 4 5],[6 7]);
A * B;

% Product of blk-vectors that results in single block (num. matrix)
A = blkmat([],[],2,[3 4 5]);
B = blkmat([],[],[3 4 5],6);
A * B;

% (Left) division by symmetric non-regular matrix
A = blkmat([],[2 3],randn(5));
B = blkmat([],[],[2 3],[4 5 6],randn(5,15));
A \ B;
% (Right) division by symmetric non-regular matrix
A = blkmat([],[],[4 5 6],[2 3],randn(15,5));
B = blkmat([],[2 3],randn(5));
A / B;

disp('All blkmat tests passed without errors')

end