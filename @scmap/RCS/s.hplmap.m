h31272
s 00010/00000/00000
d D 1.1 97/04/24 09:53:35 tad 1 0
c 
e
u
U
f b 
f e 0
t
T
I 1
function M = hplmap(M)
%HPLMAP Convert generic Schwarz-Christoffel map object to half-plane map.
%   HPLMAP(M) creates a hplmap object based on the polygon and
%   options contained in M.
%   
%   See the HPLMAP class documentation.

%   Copyright 1997 by Toby Driscoll. Last updated %G%.

M = hplmap(M.polygon,M.options);
E 1
