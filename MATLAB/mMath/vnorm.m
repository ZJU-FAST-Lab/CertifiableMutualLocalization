function m = vnorm( X )
% m = vnorm( X )
% Compute norm of column vectors in the matrix
% (vertical norm or vector norm)

m = sqrt( sum( X.^2, 1 ) );