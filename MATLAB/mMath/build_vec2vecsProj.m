function P = build_vec2vecsProj( n )
% P = build_vec2vecsProj( n )
% 
% Builds the projection matrix that fulfills the relation
%   vecs(S) = P * vec(S)
% where S is a symmetric nxn matrix,
% so P is (n(n+1)/2) x (n^2)

n_ = n*(n+1)/2;

% Create maps of indeces
idxMatVecs = avech( 1:n_ );
idxMatVec  = avec( 1:n^2 );

P = zeros(n_,n^2);
for j=1:n
  for i=j:n
    for k=1:n
      for l=1:n
        a = idxMatVecs(i,j);
        b = idxMatVec(k,l);
        if i==j && j==k && k==l
          P(a,b) = 1;
        elseif (i==k && j==l && i~=j) || (i==l && j==k && l~=j)
          P(a,b) = 1/sqrt(2);
        end
      end
    end
  end
end