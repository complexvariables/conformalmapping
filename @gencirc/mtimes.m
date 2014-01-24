function gc = mtimes(gc,z)

if isa(z,'gencirc')
  tmp = gc;
  gc = z;
  z = tmp;
end

if isinf(z)
  error('Must scale by a finite number')
end

gc.point = gc.point*z;
gc.center = gc.center*z;
gc.radius = gc.radius*abs(z);