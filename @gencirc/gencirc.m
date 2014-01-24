function gc = gencirc(varargin)

superiorto('double');

switch nargin
  case {0}
    gc.point = [];
    gc.center = [];
    gc.radius = [];
    gc.interiorpoint = [];
  case {1}
    if isa(varargin{1},'gencirc')
      gc = varargin{1};
      return
    elseif isa(varargin{1},'double') & length(varargin{1})==3
      gc.point = varargin{1}(:);
    else
      error('Must give a 3-vector of complex numbers')
    end
  otherwise
    error('Must give a 3-vector of complex numbers')
end

% Deduce center and radius
M = mobius([1 1i -1],gc.point);
A = M.matrix;
zi = feval( inv(M), Inf );  % pre-image of infinity
if abs( abs(zi) - 1 ) < 10*eps
  % This g.c. is a line
  gc.center = NaN;
  gc.radius = Inf;
  % If Inf was given, make it the last point
  if isreal(gc.point), gc.point = complex(gc.point); end
  gc.point = sort( gc.point );
else 
  % Inverse of zi must map to the center
  gc.center = M(1/conj(zi));
  gc.radius = abs( gc.point(1) - gc.center );
end

% Find a point in the interior of the curve. For a line, this is a point to
% the "left" as seen by following the given points.
if ~isinf(gc.radius)
  gc.interiorpoint = gc.center;
else
  tangent = diff(gc.point(1:2));
  gc.interiorpoint = gc.point(1) + 1i*tangent;
end

cc = closedcurve([]);
gc = class(gc,'gencirc',cc);

  