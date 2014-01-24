function eta = transpose(zeta)

eta = homog( transpose(numer(zeta)), transpose(denom(zeta)) );