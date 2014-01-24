function M = mrdivide(M1,M2)
%Divide Mobius map by a scalar, or reciprocate it.
%   1/M, for Mobius map M, swaps the numerator and denominator of M.
%   M/c, for scalar c, multiplies the denominator of M by c.

%   Copyright (c) 1998 by Toby Driscoll.
%   $Id: mrdivide.m,v 1.1 1998/07/01 20:14:22 tad Exp $

if isequal(M1,1) 
  % Exchange numerator and denominator
  A = M2.matrix;
  M = mobius(A([2 1],:));
elseif isa(M2,'double') & length(M2)==1
  A = M1.matrix;
  M = mobius( [1 0;0 M2]*A );
else
  error('Division not defined for these operands.')
end
