function e = canvec( n, k )
% e = canvec( n, k )
% Returns the canonical vector of dimension n with 1 at k-th position.

assert( k >= 1 && k <= n );

In = eye(n);
e = In(:,k);

end