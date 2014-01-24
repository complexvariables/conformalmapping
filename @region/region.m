function r = region(varargin)

% REGION(p) or REGION(p,[]) gives interior of p.
% REGION([],q) gives exterior of q.
% REGION(p,q) gives interior of p and exterior of q.

switch nargin
  case 0
    r.outerboundary = [];
    r.innerboundary = [];
    r.isin = [];
    r = class(r,'region');
    return
  case 1
    p = varargin{1};
    if isa(p,'region') 
      r = p;
    elseif isa(p,'closedcurve')
      r.outerboundary = p;
      r.innerboundary = [];
    else
      error('Single argument must be a region or a closedcurve')
    end
  otherwise
    [p,q] = deal(varargin{1:2});
    if isa(p,'closedcurve') 
      r.outerboundary = p;
    elseif isempty(p)
      r.outerboundary = [];
    else
      error('Boundary components must be closedcurves')
    end
    if isa(q,'closedcurve') 
      r.innerboundary = q;
    elseif isempty(q)
      r.innerboundary = [];
    else
      error('Boundary components must be closedcurves')
    end
end

if isempty(r.innerboundary)
  r.isin = @(z) isinside(r.outerboundary,z);
elseif isempty(r.outerboundary)
  r.isin = @(z) ~isinside(r.innerboundary,z);
else
  r.isin = @(z) isinside(r.outerboundary,z) & ~isinside(r.innerboundary,z);
end

r = class(r,'region');