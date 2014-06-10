classdef closedcurve
% CLOSEDCURVE abstract base class for simple planar Jordan curves.

% This file is a part of the CMToolbox.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.
% (Re)written by Everett Kropf, 2014,
% adapted from Toby Driscoll's code, originally 20??.

properties
  corner_                               % Placeholder.
  length_ = 1                           % Parameter length.
end

methods
  % Just use default constructor; nothing to do right now on construction.
  
  function box = boundbox(C)
    % Return bounding box for curve using evenly spaced points.
    t = length(C)*(0:199)'/200;
    box = boundbox(point(C, t));
  end
  
  function varargout = corner(C, varargin)
    % I'm honestly not quite sure about the usage of this, it will probably
    % change drastically in the future. Haven't found it being used yet in any
    % of Toby's code examples. -- EK
    % Other than being set in the polygon constructor. -- EK
    v = C.corner_;
    
    % Classes of input args?
    inclass = cellfun(@class, varargin, 'uniformoutput', false);
    
    % Did we get an index?
    j = find(strcmp('double', inclass));
    if ~isempty(j)
      v = v(varargin{j});
    end
    
    % Field selection?
    j = find(strcmp('char', inclass));
    if ~isempty(j)
      varargout{1} = cat(1, v.(varargin{j}));
    else
      if nargout == 3
        t = cat(1, v.param);
        z = cat(1, v.point);
        alpha = cat(1, v.alpha);
        varargout = {t, z, alpha};
      else
        varargout{1} = v;
      end
    end
  end
  
  function display(C)
    fprintf('\n%s =\n', inputname(1))
    disp(C)
  end
  
  function out = fill(C, varargin)
    washold = ishold;
    
    % Rely on adaptplot or overloaded plot to get curve "right".
    hold on
    h = plot(C);
    z = complex(get(h, 'xdata'), get(h, 'ydata'));
    delete(h)
    
    args = plotdef.fillargs;
    h = fill(real(z), imag(z), args{:}, varargin{:});
    
    if ~washold
      hold off
    end
    
    if nargout
      out = h;
    end
  end
  
  function n = length(C)
    % Curve parameter length.
    n = C.length_;
  end
  
  function out = plot(C, varargin)
    % Plot curve in the plane.
    washold = ishold;
    newplot
    
    h = plot_(C);
    args = plotdef.closedcurveargs;
    set(h, args{:}, varargin{:});

    if ~washold
      axis(plotbox(C, 1.1));
      set(gca, 'dataaspectratio', [1 1 1])
      hold off
    end
    
    if nargout > 0
      out = h;
    end
  end
    
  function box = plotbox(C, scale)
    % Return plot box for curve using evenly spaced points.
    %
    % See also plotbox.
    
    if nargin < 2
      scale = [];
    end
    t = length(C)*(0:199)'/200;
    box = plotbox(point(C, t), scale);
  end

  function out = rsplot(C, varargin)
    % Plot curve on the Riemann sphere.
    washold = ishold;
    
    % Draw Riemann shpere if not there.
    if isempty(findobj(gca, 'tag','CMT:RiemannSphere')) || ~washold
      [xs, ys, zs] = sphere(36);
      mesh(0.995*xs, 0.995*ys, 0.995*zs, 'edgecolor', .85*[1 1 1], ...
           'tag', 'CMT:RiemannSphere')
      hold on
    end
    
    % Draw on the sphere.
    function x = rspoint(t)
      z = point(C, t);
      [x1, x2, x3] = c2rs(z);
      x = [x1, x2, x3];
    end
    h = adaptplot(@rspoint, [0, length(C)]);
    set(h, varargin{:});
    
    if ~washold
      hold off
    end
    axis equal
    
    if nargout > 0
      out = h;
    end
  end
  
  function z = subsref(C, S)
    % Equate C(t) with point(C, t). Fallback to builtin subsref otherwise.
    if length(S) == 1 && strcmp(S.type, '()')
      if length(S.subs) == 1
        z = point(C, S.subs{1});
        return
      else
        error('Object only takes single parenthised subscript.')
      end
    end
    
    z = builtin('subsref', C, S);
  end
end

methods(Abstract=true)
  z = point(C, t)
end

methods(Hidden)
  function h = plot_(C, varargin)
    function x = xypoint(t)
      z = point(C, t);
      x = [real(z), imag(z)];
    end
    h = adaptplot(@xypoint, [0 length(C)]);
  end
end

end
