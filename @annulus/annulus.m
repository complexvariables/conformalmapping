function a = annulus(varargin)

switch nargin
  case 0
    a.center = [];
    a.innerrad = [];
    a.outerrad = [];
    a = class(a,'annulus',region);
    return
  case 2
    a.center = 0;
    [a.innerrad,a.outerrad] = deal(varargin{:});
  case 3
    [a.center,a.innerrad,a.outerrad] = deal(varargin{:});
  otherwise 
    error('Must give radii (inner first), or center and radii')
end

cin = gencirc(a.center + a.innerrad*[1 1i -1]);
cout = gencirc(a.center + a.outerrad*[1 1i -1]);
a = class(a,'annulus',between(cout,cin));
