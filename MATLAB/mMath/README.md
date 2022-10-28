# mMath toolbox
Collection of Matlab math-related (mainly algebraic) functions for improved readability and compactness.

WARNING: This repository is mainly for personal use, and some notation could change in the API.

## An overview of functions

NOTE: This list is incomplete yet

- `R = project_rotation( M )`
  Find rotation `R=U*S*V'` closest to `M=U*D*V'` minimizing the chordal distance.

- `[R,t] = project_rotation( M, fmin )`
  Find rotation `R=U*S*V'` that minimizes `fmin(R)`.

- `Z = null(A,tol)`
  The same function as Matlab but with user-definable tolerance.

- `Q = orth(A,tol)`
  The same function as Matlab but with user-definable tolerance.
