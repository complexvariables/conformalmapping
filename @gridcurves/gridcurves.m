function g = gridcurves(r,varargin)

g.region = region;
g.curve = {};

if nargin==0
  g = class(g,'gridcurves');
  return
end

if ~isa(r,'region')
  error('First argument must be a region')
else
  g.region = r;
end

% Each curve is defined by a vector of points or by a function
% parameterized over [0,1].
isfun = cellfun('isclass',varargin,'function_handle');
isdble = cellfun('isclass',varargin,'double');
if all( isfun | isdble )
  g.curve = varargin;
else
  error('Arguments after the first must be points or function handles')
end

g = class(g,'gridcurves');