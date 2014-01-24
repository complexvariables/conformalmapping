function zeta = subsref(zeta,s)

switch s.type
case '()'
  zeta = homog( subsref(numer(zeta),s), subsref(denom(zeta),s) );
otherwise
  error('This type of indexing is not supported by homog objects')
end