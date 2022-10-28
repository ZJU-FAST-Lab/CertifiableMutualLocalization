function is_it = issquare( M )
% ISSQUARE True for square matrix (2D array).

s = size(M);
if (length(s) > 2)
  is_it = false;
  return;
end
is_it = (s(1) == s(2));
end