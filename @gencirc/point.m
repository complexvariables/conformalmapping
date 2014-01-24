function z = point(gc,t)

if ~isinf(gc.radius)
  theta = angle( gc.point(1) - gc.center ) + 2*pi*t;
  z = gc.center + gc.radius*exp(1i*theta);
else
  % Use homogeneous coordinates to define a reasonable interpolant
  tangent = diff(gc.point(1:2));  % must be finite
  upper = 2*tangent*(t-1/2);
  lower = 4*t.*(1-t);
  z = double( homog(gc.point(1)) + homog(upper,lower) );
end