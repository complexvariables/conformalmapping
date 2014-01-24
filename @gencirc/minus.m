function gc = minus(gc,z)

if isinf(z)
  error('Must translate by a finite number')
end

gc.point = gc.point - z;
gc.center = gc.center - z;