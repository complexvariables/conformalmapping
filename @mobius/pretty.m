function M = pretty(M,c)
%Normalize a Moebius transformation for 'nicer' numbers.
%   PRETTY(M), where M is a moebius map, divides the coefficients of M
%   by the constant term in the demoninator, or (if that is zero) the
%   coefficient of the linear term in the denominator. Then any real or
%   imaginary parts that are close to machine precision are rounded to
%   exactly zero. 
%
%   PRETTY(M,C) instead normalizes all coefficients by C, then cleans up
%   small numbers.


%   Copyright (c) 2004, 2006 by Toby Driscoll.
%   $Id$

A = M.matrix;

% Normalization
if nargin==1
  d = A(2,2);
  if d==0, d = A(2,1); end
  A = A/d;
else
  A = c*A;
end

% Rounding
index = abs(imag(A)) < 100*eps;
A(index) = real(A(index));
index = abs(real(A)) < 100*eps;
A(index) = 1i*imag(A(index));

M = mobius(A);
