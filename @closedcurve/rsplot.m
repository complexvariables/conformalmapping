function out = rsplot(cc,varargin)

washold = ishold;
% Draw the Riemann sphere if none exists
if isempty(findobj(gca,'tag','CMT:RiemannSphere')) | ~washold
  [xs,ys,zs] = sphere(36);
  mesh(0.995*xs,0.995*ys,0.995*zs,'edgecolor',.85*[1 1 1],...
    'tag','CMT:RiemannSphere')
  hold on
end

% Do the plot
h = adaptplot( @rspoint, [0 1] );
set(h,varargin{:});

  function x = rspoint(t)
    z = point(cc,t);
    [x1,x2,x3] = c2rs(z);
    x = [ x1 x2 x3 ];
  end

% Clean up
if ~washold, hold off, end
axis equal
if nargout > 0, out = h; end
  
end