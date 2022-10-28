function out = isrotation(X)

if norm(X'*X-eye(3),'fro') > 1e-4 || norm(det(X)-1) > 1e-4
  out = false;
else
  out = true;
end