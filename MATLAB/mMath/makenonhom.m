function nonhom_x = makenonhom( x )
% nonhom_x = makenonhom( x )
% Returns the homogeneous vector (append 1 at the end of a column vector)

if size(x,2) ~= 1
  error('Use only column vectors with makenonhom');
end

if abs(x(end)) < 1e-6
  warning('Hom component is zero');
end

nonhom_x = x(1:end-1) / x(end);