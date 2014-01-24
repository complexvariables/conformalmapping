function c = mtimes(a,b)

if isfloat(a)
  a = homog(a);
elseif isfloat(b)
  b = homog(b);
end

c = homog( a.numer*b.numer, a.denom*b.denom );