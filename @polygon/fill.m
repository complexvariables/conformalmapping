function H = fill(p,varargin)
%FILL   Plot a polygon with a filled interior.
%   FILL(P) plots the boundary of P in blue and fills the interior of the
%   polygon with gray. FILL(P,PROP1,VAL1,...) passes additional arguments
%   to the built-in FILL command.
%
%   See also FILL.

% Copyright 2003 by Toby Driscoll.
% $Id: fill.m,v 1.3 2004/05/27 13:11:21 driscoll Exp $

v = vertex(p);
vf = v(~isinf(v));
if any(isinf(v))
  v = vertex(truncate(p));
end

axlim = [min(real(vf)) max(real(vf)) min(imag(vf)) max(imag(vf))];
d = max([diff(axlim(1:2)),diff(axlim(3:4))]);
if d < eps, d = 1; end
axlim(1:2) = mean(axlim(1:2)) + 0.54*[-1 1]*d;
axlim(3:4) = mean(axlim(3:4)) + 0.54*[-1 1]*d;

% Use defaults, but allow overrides and additional settings.
settings = { [0.75 0.75 0.85],'edgecolor','b','linewidth',1.5, varargin{:} };
v = v([1:end 1]);
h = fill(real(v),imag(v),settings{:});

if ~ishold
  axis equal
  axis square
  axis(axlim)
end

if nargout > 0
  H = h;
end
