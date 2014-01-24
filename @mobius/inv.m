function Minv = inv(M)
%Inverse of a Mobius transformation.

status = warning('query','MATLAB:singularMatrix');
warning('off','MATLAB:singularMatrix');
A = M.matrix;
if rcond(A) < eps
  error('Mobius transformation appears to be singular')
end
Minv = mobius(inv(A));
warning(status.state,'MATLAB:singularMatrix');
