function bool = isinside(gc,z)

if isinf(gc)
  z = (z - gc.point(1))/tangent(gc)
  bool = (imag(z) > 0);
else
  bool = abs( z-gc.center ) < gc.radius;
end
