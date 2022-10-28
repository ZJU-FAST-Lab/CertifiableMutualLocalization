function e = canvec3( k )

assert( k >= 1 && k <= 3 );

I3 = eye(3);
e = I3(:,k);

end