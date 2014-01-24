function tau = tangent(gc,z)

if isinf(gc)
  tau = diff(gc.point(1:2));
else
  tau = 1i*(z-gc.center);
end
