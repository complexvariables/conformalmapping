function r = between(p,q)

if isinside(p,point(q,0.123456789))
  r = region(p,q);
elseif isinside(q,point(p,0.123456789))
  r = region(q,p);
else
  error('Boundaries do not appear to be nested')
end