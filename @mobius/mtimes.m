function M = mtimes(M1,M2)
%Multiply Moebius transformation by a scalar.

%   Copyright (c) 1998, 2006 by Toby Driscoll.
%   $Id$

% Make the first one mobius.
if isa(M1,'double')
  tmp = M1;
  M1 = M2;
  M2 = tmp;
end

A = M1.matrix;
if isa(M2,'double') & length(M2)==1
  A(1,:) = A(1,:)*M2;
  M = mobius(A);
else
  error('Multiplication not defined for these operands.')
end
