function e = end(a,k,n)

if n==1
  e = length(a.numer);
else
  e = size(a.numer,k);
end