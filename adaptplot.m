function [out1,out2] = adaptplot(fun,tspan,varargin)
%ADAPTPLOT Adaptively plot an explicit or parametric curve.
%   ADAPTPLOT(FUN,TSPAN) adaptively selects values of t in the interval
%   TSPAN in order to produce a nice-looking plot. FUN should accept column
%   vector inputs for t and produce an array X=FUN(T), where each column of
%   X is one component of the curve. If size(X,2)==1, the plot is X versus
%   T. If size(X,2) is 2 or 3, the curve is parametrically defined by the
%   columns of X. 
%
%   ADAPTPLOT(FUN,TSPAN,'trace',DELAY) shows the iterative addition of
%   points, with an optionally specified pause between iterations.
%
%   HAN = ADAPTPLOT(FUN,TSPAN) returns a handle to the resulting line.
%
%   [T,X] = ADAPTPLOT(FUN,TSPAN) returns the computed T values and X arrays
%   without producing any graphics.
%
%   Examples:
%     adaptplot( @humps, [0,1] )
%     adaptplot( @humps, [0,1], 'trace', 0.75)
%     adaptplot( @(t) exp(-3*sin(t)+2*cos(4*t)), [2 8], 'trace', 0.75)
%     adaptplot( @(t) [cos(t) sin(t)], [-pi pi] )
%     adaptplot( @(t) [t.*cos(t), exp(t).*sin(t), t], [0 6*pi] )
%
%   See also FPLOT.

%   Copyright Toby Driscoll, 2006. All rights reserved.
%   Version 1.0, 02 June 2006

%% Preliminaries
[traceflag,delay] = optargs({'no',0.1},varargin);
dotrace = isequal(traceflag,'trace');
if dotrace, shg, end

x = feval(fun,tspan(1));
ndim = size(x,2);
if nargout < 2
  % We are expected to do some plotting. Start with a high-level plot
  % command to get the axes behavior right, then use line() to avoid
  % holding issues.
  makeplot = true;
  if ndim < 3
    allhand = plot(NaN,NaN,'k.-','erasemode','back');
  else
    allhand = plot3(NaN,NaN,NaN,'k.-','erasemode','back');
  end
  badhand = line(NaN,NaN,'color','r','marker','.','linestyle','none',...
    'erasemode','none');  
else
  makeplot = false;
end

fun = fcnchk(fun);  % this is in FPLOT

% If the axis limits are frozen, use them to determine the stopping
% criterion. Otherwise, it will be based on the size of the graph itself.
if strcmp(get(gca,'xlimmode'),'manual') & ishold
  axlim = axis;
  axlim = reshape( axlim, [2 length(axlim)/2] );
  diam = max( diff(axlim,[],1) );
  axfixed = true;
else
  diam = 0;
  axfixed = false;
end

maxdepth = 12;  % number of iterations of refinement
n = 14;         % initial number of points (prime # of intervals)
tol = 2e-3;     % acceptable error relative to diam

%% Begin refinement
t = linspace(tspan(1),tspan(2),n)';
x = feval(fun,t);
refine = true;
depth = 1;
tnew = [];  xnew = zeros(0,ndim);

while any(refine) & (depth < maxdepth)
  
  % Update trace plot
  if makeplot & dotrace
    updateplot(allhand,t,x,ndim)
    updateplot(badhand,tnew,xnew,ndim)
    drawnow, pause(delay)
  end
  
  % Perform forward/backward linear extrapolation from each pair of
  % neighboring points
  dx = diff(x,[],1);  dt = diff(t);
  xf = x(1:n-2,:) + dx(1:n-2,:) .* ...
    repmat( (t(3:n)-t(1:n-2))./dt(1:n-2), [1 ndim] );
  xb = x(3:n,:) + dx(2:n-1,:) .* ...
    repmat( (t(1:n-2)-t(3:n))./dt(2:n-1), [1 ndim] );
  
  % Errors in forward and backward differences
  ef = ptdist( xf, x(3:n,:) );
  eb = ptdist( xb, x(1:n-2,:) );

  % Make refinement decisions
  refine = false(n,1);
  refine(1:n-2) = (eb > tol*diam);
  refine(3:n) = refine(3:n) | (ef > tol*diam);
  
  % Offending points add new t values 1/3 of the way toward either side
  ref1 = find( refine(1:n-1) );
  ref2 = find( refine(2:n) );
  tnew = [ t(ref1) + dt(ref1)/3; t(1+ref2) - dt(ref2)/3 ];
  xnew = feval(fun,tnew);
  
  % Sort in the new entries
  [t,idx] = sort( [t;tnew] );
  x = [x;xnew];  x = x(idx,:);

  % Update parameters
  n = length(t);
  depth = depth+1;
  if ~axfixed
    diam = abs( max( max(x,[],1) - min(x,[],1) ) );
  end
end

%% Wrap up
if makeplot
  set(allhand,'erasemode','normal','marker','none','color','b')
  updateplot(allhand,t,x,ndim)
  set(badhand,'erasemode','normal')
  delete(badhand)
end

if nargout==2
  out1 = t;  out2 = x;
elseif nargout==1
  out1 = allhand;
end

end

%% Subfunctions

function d = ptdist(x,y)
d = sqrt( sum( abs(x-y).^2, 2 ) );
end

function updateplot(hand,t,x,ndim)
if ndim==1
  set(hand,'xdata',t,'ydata',x)
elseif ndim==2
  set(hand,'xdata',x(:,1),'ydata',x(:,2))
else
  set(hand,'xdata',x(:,1),'ydata',x(:,2),'zdata',x(:,3))
end
end

function varargout = optargs(default,arg)
% Use arg to override defaults
varargout = default;
idx = find( ~cellfun('isempty',arg) );
varargout(idx) = arg(idx);
end

    