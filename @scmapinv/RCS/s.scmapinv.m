h37729
s 00012/00000/00000
d D 1.1 97/07/22 14:41:35 tad 1 0
c 
e
u
U
f b 
f e 0
t
T
I 1
function Mi = scmapinv(M)
%SCMAPINV Inverse of a Schwarz-Christoffel map.
%   SCMAPINV(M) returns a dummy object that represents the inverse of
%   the SC map M. The only things possible with this object are to EVAL
%   it or INV it back to the SC map.
%   
%   See also SCMAPINV/EVAL, SCMAPINV/SUBSREF, SCMAPINV/INV.

%   Copyright 1997 by Toby Driscoll. Last updated %G%.

Mi.themap = M;
Mi = class(Mi,'scmapinv');
E 1
