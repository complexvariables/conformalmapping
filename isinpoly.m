function in = isinpoly(z,vertex)

in = inpolygon(real(z),imag(z),real(vertex),imag(vertex));
