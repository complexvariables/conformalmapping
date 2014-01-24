function c = plus(a,b)

if isfloat(a)
  a = homog(a);
end
if isfloat(b)
  b = homog(b);
end

c = homog( a.numer*b.denom + a.denom*b.numer, a.denom*b.denom );