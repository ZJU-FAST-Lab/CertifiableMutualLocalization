function B = subsref(A,S)

if length(S)==1 % A(i,j) 
  switch S.type
   case '()', 
    [rows, cols] = extract_indices(A, S.subs);
    if max(rows) > nrows(A) | max(cols) > ncols(A)
      error(['index ' num2str(rows) ', ' num2str(cols) ' is out of bounds']);
    end
    B = A.M(block(rows, rowsizes(A)), block(cols, colsizes(A)));
   otherwise,
    error(['Unrecognized subscript operator ' S.type]);
  end
else
  error(['Unrecognized subscript operator ' S]);
end
