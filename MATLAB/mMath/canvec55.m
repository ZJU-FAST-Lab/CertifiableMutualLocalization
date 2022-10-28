function E = canvec55( i,j )

dim = 5;

assert( i >= 1 && i <= 5 );
assert( j >= 1 && j <= 5 );

E = zeros(dim,dim);
E(i,j) = 1;

end