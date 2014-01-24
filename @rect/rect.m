function r = rect(varargin)

switch nargin
  case 0
    r.corner = [];
    r = class(r,'rect',region);
    return
  case 1
    boundary = polygon(varargin{1});
    r.corner = vertex(boundary);
  case 2
    [width,height] = deal(varargin{:});
    boundary = polygon( [0; width; width+1i*height; 1i*height] );
    r.corner = vertex(boundary);
  otherwise 
    error('Must give a boundary tangent, or point and tangent')
end

r = class(r,'rect',interior(boundary));
