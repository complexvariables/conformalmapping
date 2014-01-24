function c = curve(g,n)

if nargin==1
  c = g.curve;
else
  c = g.curve{n};
end
