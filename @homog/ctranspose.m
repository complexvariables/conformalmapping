function eta = ctranspose(zeta)

eta = homog( ctranspose(numer(zeta)), ctranspose(denom(zeta)) );