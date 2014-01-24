function Mi = scmapinv(M)
%SCMAPINV Inverse of a Schwarz-Christoffel map.
%   SCMAPINV(M) returns a dummy object that represents the inverse of
%   the SC map M. The only things possible with this object are to EVAL
%   it or INV it back to the SC map.
%   
%   See also SCMAPINV/EVAL, SCMAPINV/SUBSREF, SCMAPINV/INV.

%   Copyright 1998 by Toby Driscoll.
%   $Id: scmapinv.m,v 2.2 2001/07/20 13:53:44 driscoll Exp $

Mi.originalmap = [];

switch nargin
  case 0
    Mi = class(Mi,'scmapinv',conformalmap);
    return
  case 1
    if isa(M,'scmapinv')
      Mi = M;  return
    elseif isa(M,'scmap')
      Mi.originalmap = M;
      Mi = class(Mi,'scmapinv',conformalmap([],image(M),domain(M)));
    end
end
