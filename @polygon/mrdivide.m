function r = mrdivide(p,q)
%   Divide a polygon by a scalar.

%   Copyright 2003 by Toby Driscoll.
%   $Id: mrdivide.m,v 1.1 2003/03/03 16:17:37 driscoll Exp $

if ~isa(q,'double') || length(q)>1
  error('Function ''/'' defined only for a scalar double.')
end

r = p;
r.vertex = r.vertex/q;
