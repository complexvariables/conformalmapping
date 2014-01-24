function [z,idx] = point(p,t)

%   Copyright 1998-2002 by Toby Driscoll. 
%   $Id: linspace.m,v 2.4 2002/09/05 15:40:11 driscoll Exp $

n = length(p);
zc = p.vertex;
zh = p.hvertex;
zhind = p.hindex;
z = repmat(NaN,size(t));
t = t(:);

% Loop by side number
for k = 1:n
  sidek = find( t>=(k-1)/n & t<=k/n );
  if isempty(sidek), continue, end
  tk = n*t(sidek) - (k-1);
  k1 = mod(k,n) + 1;
  if isinf(zc(k))
    tangent = numer( zh(zhind(k)+1) );
    z(sidek) = double( zc(k1) + homog( (1-tk)*tangent, tk ) );
  elseif isinf(zc(k1))
    tangent = numer( zh(zhind(k1)) );
    z(sidek) = double( zc(k) + homog( tk*tangent, 1-tk ) );
  else
    z(sidek) = interp1( [0 1], zc([k k1]), tk, 'linear' );
  end
end



