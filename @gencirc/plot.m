function out = plot(gc,varargin)

newplot

if isinf(gc)
  % For a line, we will use a polygon. This employs the truncation
  % mechanism that gives us something usable with plotting regions.
  h = plot( polygon( [gc.point(1) infvertex(tangent(gc),-tangent(gc))] ));
  
%   % For a line, we can plot 2 points
%   if isequal( get(gca,'xlimmode'), 'auto' ) 
%     % These are quite arbitrary
%     z = gc.point(1) + tangent(gc)*[-2 2];
%   else
%     % We need to show the part of the line inside the current axes box. We
%     % do it by intersecting with the circumscribing circle.
%     axlim = axis;
%     circbox = gencirc( axlim([1 2 2]) + 1i*axlim([1 1 2]) );
%     z = intersect(gc,circbox);
%   end
%   h = plot(real(z),imag(z));
else
  % Plot a genuine circle
  h = adaptplot( @realpoint, [0 1] );
  % Plot a slightly larger invisible box around to make the axis limits
  % look nicer.
  xl = real(gc.center) + gc.radius*[-1.02 1.02];
  yl = imag(gc.center) + gc.radius*[-1.02 1.02];
  line( xl([1 2 2 1]), yl([1 1 2 2]), 'visible','off',...
    'handlevisibility','off' );
end

% Clean up
set(h,'linewidth',1.25,varargin{:});
axis image
if nargout > 0, out = h; end

  function x = realpoint(t)
    z = point(gc,t);
    x = [ real(z) imag(z) ];
  end

end