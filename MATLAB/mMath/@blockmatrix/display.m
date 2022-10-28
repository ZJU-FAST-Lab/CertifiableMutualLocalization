function display(A)

if A.row_regular & A.col_regular
  s = sprintf('%dx%d matrix of %dx%d blocks', nrows(A), ncols(A), rowsize(A), colsize(A));
elseif A.row_regular & ~A.col_regular
  s = sprintf('%dx- matrix of %d x [%s] blocks', nrows(A), rowsize(A), num2str(colsizes(A)));
elseif ~A.row_regular & A.col_regular
  s = sprintf('-x%d matrix of [%s] x %d blocks', ncols(A), num2str(rowsizes(A)), colsize(A));
else
  s = sprintf('-x- matrix of [%s] x [%s] blocks', num2str(rowsizes(A)), num2str(colsizes(A)));
end
disp(' ');
disp([' ' inputname(1) ' = ' s ]);
disp(' ');
disp(A.M)
