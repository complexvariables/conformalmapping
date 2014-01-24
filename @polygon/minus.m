function r = minus(p,q)
%   Translate a polygon, or subtract the vertices of two polygons.

%   Copyright 1999-2003 by Toby Driscoll.
%   $Id: minus.m,v 2.2 2003/03/03 16:28:22 driscoll Exp $

r = plus(p,-q);
