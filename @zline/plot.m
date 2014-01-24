function out = plot(this,varargin)

newplot

% We will use a polygon. This employs the truncation
% mechanism that gives us something usable for plotting the interior.

tau = tangent(this);
h = plot( polygon( [this.base infvertex(tau,-tau)] ));

% Clean up
set(h,varargin{:});
if nargout > 0, out = h; end

end