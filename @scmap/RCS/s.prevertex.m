h27117
s 00010/00000/00000
d D 1.1 97/09/18 13:42:05 tad 1 0
c 
e
u
U
f b 
f e 0
t
T
I 1
function z = prevertex(M)
%PREVERTEX Extract a vector of the prevertices of an S-C map.

tmp = parameters(M);
if strmatch('prevertex',fields(tmp))
  z = tmp.prevertex;
else
  msg = sprintf('Prevertices not defined for map of class %s\n',class(M));
  error(msg)
end
E 1
