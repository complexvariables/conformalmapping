function z = intersect(gc1,gc2)

% Map the first circle to the real axis
M = mobius(point(gc1,[1/3 2/3 1]),[-1 1 Inf]);
gc = M(gc2);
if isinf(gc)
  % Intersect real axis with a line
  tau = tangent(gc);  p = gc.point(1);
  if abs(imag(tau)) > 100*eps 
    t = -imag(p)/imag(tau);
    z = real(p) + t*real(tau);
  else 
    warning('CMToolbox:illConditionedIntersection',...
      ['Circles are close to tangency.\nIntersection problem'...
      ' is not well conditioned.'])
    z = [];
  end
  z = [z Inf];
else
  % Intersect real axis with a circle
  rat = -imag(gc.center)/gc.radius;
  if abs( abs(rat) - 1 ) < 100*eps    
    warning('CMToolbox:illConditionedIntersection',...
      ['Circles are close to tangency.\nIntersection problem'...
      ' is not well conditioned.'])
  end
  theta = asin(rat);    % find one intersection
  theta = theta(isreal(theta));             % may not have one
  theta = unique( [theta pi-theta] );       % may have a second
  z = real( gc.center + gc.radius*exp(1i*theta) );
end
z = feval(inv(M),z);
