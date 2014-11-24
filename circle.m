classdef circle < closedcurve
% CIRCLE is a generalized circle class.
%
% C = circle(center, radius)
%   Creates a circle with given center and radius.
%
% C = circle([z1, z2, z3])
%   Creates a generalized circle passing through the three given points.

% This file is a part of the CMToolbox.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.
% (Re)written by Everett Kropf, 2014,
% adapted from code by Toby Driscoll, originally 20??.

properties
  points
  circCenter
  circRadius
  interiorPoint
end

properties(Hidden,Constant)
    ei2pi = exp(2i*pi*(0:199)'/199)
end

methods
  function gc = circle(varargin)
    if ~nargin
      return
    end
    
    if isa(varargin{1}, 'circle')
      gc = varargin{1};
      return
    end
    
    badargs = true;
    switch nargin
      case 1
        z3 = varargin{1};
        if isa(z3, 'double') && numel(z3) == 3
          badargs = false;
          
          % Deduce center and radius.
          % Use standardmap directly to avoid vicious circularity, since mobius
          % constructs a circle when given two 3-tuples.
          M = mobius(mobius.standardmap(z3)\mobius.standardmap([1, 1i, -1]));
          zi = pole(M);
          if abs(abs(zi) - 1) < 10*eps
            % Infinite radius case.
            center = nan;
            radius = inf;
            % If Inf was given, make it the last point.
            if isreal(z3)
              z3 = complex(z3);
            end
            z3 = sort(z3);
          else
            % Inverse of zi maps to center.
            center = M(1/conj(zi));
            radius = abs(z3(1) - center);
          end
          
          % Find a point in the interior of the curve. For a line, this is a
          % point to the "left" as seen by following the given points.
          if ~isinf(radius)
            interiorpt = center;
          else
            tangent = diff(z3(1:2));
            interiorpt = z3(1) + 1i*tangent;
          end
        end
        
      case 2
        [center, radius] = varargin{:};
        centercond = isa(center, 'double') && numel(center) == 1 ...
          && ~(isnan(center) || isinf(center));
        radiuscond = isa(radius, 'double') && numel(radius) == 1 ...
          && ~(isnan(radius) || isinf(radius) || radius < 0);
        if (centercond && radiuscond)
          badargs = false;
          
          if radius > 0
            z3 = center + radius*exp(1i*pi*[0, 0.5, 1]);
          else
            z3 = []; % Degenerate circle.
          end
          interiorpt = center;
        end
    end
    if badargs
      error('Circle takes a vector of 3 points or a center and radius.')
    end
    
    gc.points = z3;
    gc.circCenter = center;
    gc.circRadius = radius;
    gc.interiorPoint = interiorpt;
  end
  
  function gc = apply(gc, m)
      if ~isa(m, 'mobius')
          error('CMT:NotDefined', ...
              'Expected a mobius transformation.')
      end
      
      gc = circle(m(gc.points));
  end
  
  function z = center(gc)
    z = gc.circCenter;
  end
  
  function disp(gc)
    if isinf(gc)
      fprintf('circle (generalized) as a line,\n')
    else
      fprintf('circle with center %s and radius %s,\n', ...
              num2str(gc.circCenter), num2str(gc.circRadius))
    end
    if isempty(gc.points)
      fprintf('\n(degenerate circle)\n\n')
    else
      fprintf('\npassing through points:\n\n')
      disp(gc.points(:))
    end
  end
  
  function d = dist(gc, z)
    % Distance between point and circle.
    if ~isinf(gc)
      v = z - gc.circCenter;
      d = abs(abs(v) - gc.circRadius);
    else
      v = z - gc.points(1);
      s = sign(1i*diff(gc.points(1:2)));
      d = abs(real(v)*real(s) + imag(v)*imag(s));
    end
  end
  
  function out = fill(gc, varargin)
      args = cmtplot.fillargs;
      z = gc.center + gc.radius*gc.ei2pi;
      h = fill(real(z), imag(z), args{:}, varargin{:});      
      if nargout
          out = h;
      end
  end
  
  function z = intersect(gc1, gc2)
    % Calculate circle intersections.
    
    % Map first circle to the real axis.
    M = mobius(point(gc1, [1/3, 2/3, 1]), [-1, 1, Inf]);
    gc = M(gc2);
    if isinf(gc)
      % Intersect real axis with a line.
      tau = tangent(gc);
      p = gc.points(1);
      if abs(imag(tau)) > 100*eps
        t = -imag(p)/imag(tau);
        z = real(p) + t*real(tau);
      else
        warning(['Circles are close to tangency.\nIntersection', ...
                 ' problem is not well conditioned.'])
        z = [];
      end
      z = [z, inf];
    else
      % Intersect real axis with a circle.
      rat = -imag(gc.circCenter)/gc.circRadius;
      if abs(abs(rat) - 1) < 100*eps
        warning(['Circles are close to tangency.\nIntersection', ...
                 ' problem is not well conditioned.'])
      end
      theta = asin(rat);                    % find one intersection
      theta = theta(isreal(theta));         % may not have one
      theta = unique([theta, pi - theta]);  % may have a second
      z = real(gc.circCenter + gc.circRadius*exp(1i*theta));
    end
    z = feval(inv(M), z);
  end
  
  function tf = isinf(gc)
    tf = isinf(gc.circRadius);
  end
  
  function tf = isinside(gc, z)
    if isinf(gc)
      z = (z - gc.points(1))/tangent(gc, z);
      tf = imag(z) > 0;
    else
      tf = abs(z - gc.circCenter) < gc.circRadius;
    end
  end
  
  function gc = minus(gc, z)
      if isa(z, 'circle')
          gc = plus(-gc, z);
      else
          gc = plus(gc, -z);
      end
  end
  
  function gc = mtimes(gc, z)
    if isa(z, 'circle')
      [z, gc] = deal(gc, z);
    end
    
    if isinf(z)
      error('Must scale by a finite number.')
    end
    
    gc.points = gc.points*z;
    gc.circCenter = gc.circCenter*z;
    gc.circRadius = gc.circRadius*abs(z);
  end
  
  function gc = plus(gc, z)
    if isa(z, 'circle')
      [z, gc] = deal(gc, z); 
    end
    if isinf(z)
      error('Must translate by a finite number.')
    end
    
    gc.points = gc.points + z;
    gc.circCenter = gc.circCenter + z;
  end
  
  function z = point(gc, t)
    if ~isinf(gc.circRadius)
      theta = angle(gc.points(1) - gc.circCenter) + 2*pi*t;
      z = gc.circCenter + gc.circRadius*exp(1i*theta);
    else
      % Use homogeneous coordinates to define a reasonable interpolant.
      tangent = diff(gc.points(1:2)); % must be finite
      upper = 2*tangent*(t - 1/2);
      lower = 4*t.*(1 - t);
      z = double(homog(gc.points(1)) + homog(upper, lower));
    end
  end
    
  function r = radius(gc)
    r = gc.circRadius;
  end
  
  function zt = tangent(gc, t)
    if isinf(gc)
      zt = diff(gc.points(1:2));
    else
      zt = 1i*(point(gc, t) - gc.circCenter);
    end
  end
  
  function gc = uminus(gc)
    gc.points = -gc.points;
    gc.circCenter = -gc.circCenter;
  end
end

methods(Hidden)
  function h = plotCurve(gc)
    if isinf(gc)
      % For a line, we will use a polygon. This employs the truncation
      % mechanism that gives us something usable with plotting regions.
      h = plot(polygon([gc.points(1), infvertex(tangent(gc), -tangent(gc))]));
    else
      % Circle!!
      z = gc.center + gc.radius*gc.ei2pi;
      h = plot(real(z), imag(z));
    end
  end  
end

end
