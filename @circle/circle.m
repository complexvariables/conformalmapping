function this = circle(varargin)

this.center = [];
this.radius = [];

switch nargin
  case 1
    if isa(varargin{1},'circle')  
      this = varargin{1}; return  % self-return
    end
    r = varargin{1};
    % Radius, or three points in a vector?
    if isa(r,'double') & length(r)==1
      this.center = 0;
      this.radius = r;
    elseif isa(r,'double') & length(r)==3
      % Check for collinearity.
      z1 = r(2)-r(1);  z2 = r(3)-r(1);
      area = abs( imag( z1'*z2 ) );
      if area < eps*max(abs(z1),abs(z2))
        warning('Circle points are numerically collinear')
      end
      
      % Find the center through mobius maps and inverse points.
      F = mobius([1 1i -1],r);
      zp = pole(F);
      this.center = F( 1/conj(zp) );
      this.radius = mean( abs( r - this.center ) );
    else
      error('Expected a radius or a vector of three points')
    end
  case 2
    % Look for a center and radius
    [zc,r] = deal(varargin{:});
    if isa(zc,'double') & isa(r,'double') & length(zc)==1 & length(r)==1
      this.center = zc;
      this.radius = r;
    else
      error('Expected a center and radius')
    end
  otherwise
    error('Too many input arguments')
end

this = class(this,'circle',closedcurve);

      