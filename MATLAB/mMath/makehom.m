function hom_x = makehom( x )
% hom_x = makehom( x )
% Returns the homogeneous vector (append 1 at the end of a column vector)

m = size(x,2);
hom_x = [x; ones(1,m)];