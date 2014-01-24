function out = plot(this,varargin)

newplot
washold = ishold;

if isempty(get(gca,'children'))
  % Find a clip box for the data.
  axbox = clippingbox( point(this,0:0.05:1), 1.1 );
else
  axbox = axis;
end
axis equal, axis(axbox)
box on
hold on

% Do the plot
h = adaptplot( @xypoint, [0 1] );
set(h,'linewidth',0.75,varargin{:});

  function x = xypoint(t)
    z = point(this,t);
    x = [ real(z) imag(z) ];
  end

% Clean up
if ~washold, hold off, end
axis equal
if nargout > 0, out = h; end
  
end