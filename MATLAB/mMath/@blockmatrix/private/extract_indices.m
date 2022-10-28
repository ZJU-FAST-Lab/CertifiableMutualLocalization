function [rows, cols] = extract_indices(A, subs)

if length(subs)==1
  if ncols(A)==1 % col vector
    rows = subs{1};
    cols = 1;
  elseif nrows(A) == 1 % row vector
    cols = subs{1};
    rows = 1;
  else
    error('must specify two indices');
  end
else
  rows = subs{1};
  cols = subs{2};
end

 % just using the test "rows == ':'" produces strange results if rows is 58, since char(58)==':'
if isstr(rows) & rows == ':'
  rows = 1:nrows(A);
end

if isstr(cols) & cols == ':'
  cols = 1:ncols(A);
end


