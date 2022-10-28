function R = LangevinId(info)
% R = LangevinId(info)
%   Generate values from a non-isotropic Langevin distribution around the
%   mean identity matrix:
%     pdf: ... * exp( trace( info * R ) )
%   If info is a nxnxN matrix, generate nxnxN output with N samples.
%  
% See also rnd.gauss.
% 
% Jesus Briales, University of Malaga, June 17, 2016.

N = size(info,3);
assert(issymmetric(info(:,:,1)),'Info must be a square symmetric matrix');
n = size(info,1);
% Compute offsets (upper bound on tr(info*Z)):
offset = multitrace(info);

% Get random samples using a naive rejection scheme
R = zeros(n, n, N);
todo = 1:N;
numTodo = length(todo);
while numTodo > 0
  % TODO: The rotations could be sampled from a tangent-Gaussian
  % then the rejection test be modified by the proposal PDF
  % This should result in increased efficiency for appropiate proposal PDF
  
  % Sample as many random rotations as TODO yet, uniformly
  R(:, :, todo) = randrot(n, numTodo);
  % Compare the function prop. to the PDF (<=1) to a random seed (0-1)
  % If the PDF value for an R is very high, it is likely that rand < PDF
  todoYet = rand(numTodo, 1) >= ...
    exp( multitrace(multiprod(info(:,:,todo),R(:,:,todo)))-offset(todo) );
  todo = todo(todoYet);
  numTodo = length(todo);
end

end