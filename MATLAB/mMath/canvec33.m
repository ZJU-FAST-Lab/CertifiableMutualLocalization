function E = canvec33( i,j )

dim = 3;

assert( i >= 1 && i <= dim );
assert( j >= 1 && j <= dim );

E = zeros(dim,dim);
E(i,j) = 1;

end