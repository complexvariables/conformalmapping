function val = isinside(p,z)

v = p.vertex;
val = inpolygon( real(z),imag(z), real(v),imag(v) );
