function q = uminus(p)
%   Negate the vertices of a polygon.
%   This may have surprising consequences if p is unbounded.

%   Copyright 2003 by Toby Driscoll (driscoll@math.udel.edu).
%   $Id: uminus.m,v 1.1 2003/03/03 16:28:04 driscoll Exp $

q = polygon( -vertex(p), angle(p) );

