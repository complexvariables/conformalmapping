function w = subsref(gc,S)
%SUBSREF Evaluate or compose maps.
%   M(Z), where M is a Mobius transformation and Z is a vector of
%   points, returns the image of the points in Z. This just a synonym
%   for FEVAL(M,Z). 
%   
%   See also MOBIUS/FEVAL, MOBIUS.

%   Copyright (c) 1998, 2006 by Toby Driscoll.
%   $Id$

if length(S) == 1 & strcmp(S.type,'()')
  error('Syntax not supported')
elseif length(S) == 1 & strcmp(S.type,'.')
  if strcmp(S.subs,'point')
    w = gc.point;
  else
    error('Field not available')
  end
else
  error('Only syntax for GENCIRC is a single parenthesized subscript.')
end
