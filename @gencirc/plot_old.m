function out = plot(gc,varargin)

r = gc.radius;
c = gc.center;
if ~isinf(r)
  z = c + r*exp(1i*linspace(-pi,pi,120));
  h = plot(z,varargin{:});
%  h = rectangle('position',[real(c)-r imag(c)-r 2*r 2*r],...
%    'curvature',[1 1],'edgecolor',[0 0 1]);
else
  ax = axis;
  diam = max( ax(2)-ax(1), ax(4)-ax(3) );
  % Note: first 2 points are guaranteed to be finite
  tangent = diam*sign( diff(gc.point(1:2)) );
  h = plot( gc.point(1)+tangent*[-1;1+eps*i], varargin{:} );
end

axis equal
if nargout > 0, out = h; end
  