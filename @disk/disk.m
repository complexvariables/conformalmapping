function r = disk(varargin)

switch nargin
  case 0
    r.center = [];
    r.radius = [];
    r = class(r,'disk',region);
    return
  case 1
    r.center = 0;
    r.radius = varargin{1};
  case 2
    [r.center,r.radius] = deal(varargin{:});
  otherwise 
    error('Must give a radius, or center and radius')
end

if ~isempty(r.center)
  circ = gencirc(r.center + r.radius*[1 1i -1]);
end

r = class(r,'disk',interior(circ));
