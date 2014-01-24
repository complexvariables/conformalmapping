function hp = halfplane(varargin)

switch nargin
  case 0
    hp.point = [];
    hp.tangent = [];
    hp = class(hp,'halfplane',region);
    return
  case 1
    hp.point = 0;
    hp.tangent = varargin{1};
  case 2
    [hp.point,hp.tangent] = deal(varargin{:});
  otherwise 
    error('Must give a boundary tangent, or point and tangent')
end

if ~isempty(hp.point)
  boundary = gencirc(hp.point + hp.tangent*[0 1 Inf]);
end

hp = class(hp,'halfplane',interior(boundary));
