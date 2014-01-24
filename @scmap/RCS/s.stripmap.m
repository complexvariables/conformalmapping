h32087
s 00010/00000/00000
d D 1.1 97/04/24 09:55:16 tad 1 0
c 
e
u
U
f b 
f e 0
t
T
I 1
function M = stripmap(M)
%STRIPMAP Convert generic Schwarz-Christoffel map object to strip map.
%   STRIPMAP(M) creates a stripmap object based on the polygon and
%   options contained in M.
%   
%   See the STRIPMAP class documentation.

%   Copyright 1997 by Toby Driscoll. Last updated %G%.

M = stripmap(M.polygon,M.options);
E 1
