function X = avecs(x)
% X = avecs(x)
% Inverse of symmetric (half) vectorization for symmetric matrices
% 
% Scaled half vectorization so that tr(ST)=vecs(S)'*vecs(T)
% Taken from "Schäcke, K. (2013). On the kronecker product."

k = numel(x);
s = -0.5 + sqrt(0.25+2*k);
X = zeros(s);

% Insert vector to lower triangular block
X( tril(true(s)) ) = x;

% Scale the off-diagonal lower triangular block with 1/sqrt(2)
X( ~triu(true(s)) ) = X( ~triu(true(s)) ) * 2 / sqrt(2);

% Symmetrize matrix (x2 factor above to sum both terms
X = symmetrize( X );

end