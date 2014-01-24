function z = intersect(zl1,zl2)

% Set up 2x2 linear system for the intersection
A = [ real(zl1.tangent) real(zl2.tangent); 
  imag(zl1.tangent) imag(zl2.tangent)];
dz = zl2.point - zl1.point;
b = [ real(dz); imag(dz) ];

if rcond(A) < eps
  warning('Lines are parallel to working precision')
  Q = null(A);
  if norm(Q'*b) < eps*norm(b)
    % Lines are essentially the same, so return any point
    z = zl1.point;
  else
    z = [];
  end
else
  p = A\b;
  z = zl1.point + p(1)*zl1.tangent;
end
