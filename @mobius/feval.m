function f = feval(this,z)
%Evaluate Mobius transformation.

%   Copyright (c) 2006 by Toby Driscoll.

switch(class(z))
  case {'circle','zline'}
    zp = pole(this);
    if dist(z,zp) < 10*eps*abs(zp)
      % Result appears to be a line.
      zp = feval( this, point(z,[1 3]/4) );
      f = zline(zp);
    else
      % Find new circle using three points.
      zp = feval( this, point(z,[1 2 3]/4) );
      f = circle(zp);
    end
  case 'double'
    f = NaN*zeros(size(z));
    
    % Convert inputs to homogeneous coordinates, and reshape
    z = homog(z);
    Z = [ numer(z(:)).'; denom(z(:)).' ];
    
    % Apply map
    W = this.matrix*Z;
    
    % Convert to complex without DBZ warnings
    f(:) = double( homog(W(1,:),W(2,:)) );
  otherwise
    error('mobius maps can be applied to floats or circles/zlines only')
end