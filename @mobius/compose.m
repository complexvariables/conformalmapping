function M = compose(varargin)
%Compose Moebius transformations.
%   COMPOSE(M1,M2,...) returns a single Moebius transformation that
%   composes all the given ones. The last given transformation is the first
%   to be applied, and vice versa.
%
%   See also MOEBIUS.

%   Copyright (c) 2006 by Toby Driscoll.
%   $Id$

A = varargin{end}.matrix;
for n = nargin-1:-1:1
  A = varargin{n}.matrix*A;
end
M = mobius(A);
