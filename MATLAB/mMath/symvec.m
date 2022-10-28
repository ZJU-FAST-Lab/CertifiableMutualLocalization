function x = symvec(var,n,varargin)
% x = SYMVEC(var,n)
% Create a symbolic vector of basename var with n elements
% 
% SYMVEC(___, set) sets the assumption that the variables belong to a set.
% Here, set can be 'real', 'positive', 'integer', or 'rational'.
% 
% See also SYM, SYMMAT

x = sym([var,'_%d'],[n 1],varargin{:});