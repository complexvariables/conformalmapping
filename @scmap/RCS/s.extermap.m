h32343
s 00010/00000/00000
d D 1.1 97/04/24 09:53:23 tad 1 0
c 
e
u
U
f b 
f e 0
t
T
I 1
function M = extermap(M)
%EXTERMAP Convert generic Schwarz-Christoffel map object to exterior map.
%   EXTERMAP(M) creates a extermap object based on the polygon and
%   options contained in M.
%   
%   See the EXTERMAP class documentation.

%   Copyright 1997 by Toby Driscoll. Last updated %G%.

M = extermap(M.polygon,M.options);
E 1
