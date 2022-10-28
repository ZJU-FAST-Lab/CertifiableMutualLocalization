function A = subsasgn(A,S,B)

if length(S)==1 % A(i,j)
  switch S.type
   case '()', 
    [rows, cols] = extract_indices(A, S.subs);
    
    % If we reference a non-existent cell, expand the array if regular
    r = max(rows);
    if r > nrows(A)
      if A.row_regular
	A.rsizes = repmat(A.rsizes(1), 1, r);
      else
	error('can''t expand row irregular blockmatrix');
      end
    end
    c = max(cols);
    if c > ncols(A)
      if A.col_regular
	A.csizes = repmat(A.csizes(1), 1, c);
      else
	error('can''t expand column irregular blockmatrix');
      end
    end
    
    scalar_rows = block(rows, A.rsizes);
    scalar_cols = block(cols, A.csizes);
    A.M(scalar_rows, scalar_cols) = B;
   otherwise,
    error(['Unrecognized subscript operator ' S.type]);
  end
else 
  error(['Unrecognized subscript operator ' S]);
end
