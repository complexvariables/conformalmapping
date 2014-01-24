h41514
s 00001/00001/00010
d D 1.2 97/07/22 14:52:38 tad 2 1
c Typo
c 
e
s 00011/00000/00000
d D 1.1 97/07/22 14:40:53 tad 1 0
c 
e
u
U
f b 
f e 0
t
T
I 1
function zp = eval(Mi,wp)
D 2
%SUBSREF Evaluate inverse SC map.
E 2
I 2
%EVAL Evaluate inverse SC map.
E 2
%   EVAL(MI,WP), where MI is an SCMAPINV object and WP is a vector of points
%   in the polygon of the map, returns the inverse image of WP under the
%   map (the forward image under MI).
%   
%   See also SCMAPINV, SCMAPINV/SUBSREF.

%   Copyright 1997 by Toby Driscoll. Last updated %G%.

zp = evalinv(Mi.themap,wp);
E 1
