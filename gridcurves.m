classdef gridcurves
% GRIDCURVES class holds grid curves for regions.
%
% grid = gridcurves(curves)
% The given cell array curves is stored in the object. This is mainly to
% facilitate the use of plot(grid) to give a standard look to grid plots.
%
% Stores grid curves as entries in cell array. Consider
%    gd = gridcurves;
% We overload gd(n) and gd(n,m) to retrieve those cell array entries. Why
% not just use the '{}' syntax? Wouldn't it be clearer we're using cell arrays?

% This file is a part of the CMToolbox.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.
% (Re)written by Everett Kropf, 2014,
% adapted from an idea by Toby Driscoll, 20??.

properties
  curves_
end

methods
  function gc = gridcurves(curves)
    if ~nargin
      return
    end
    
    if nargin > 1 || ~isa(curves, 'cell')
      error('CMT:InvalidArgument', ...
            'Expected a cell array of individual grid curves.')
    end
    
    gc.curves_ = curves;
  end
  
  function gc = conj(gc)
    % Complex conjugation.
    for k = 1:numel(gc.curves_)
      gc.curves_{k} = conj(gc.curves_{k});
    end
  end
  
  function disp(gd)
    fprintf('gridcurves object:\n\n')
    fprintf('  with %d gridlines.\n\n', numel(gd.curves_))
  end
  
  function gc = minus(gc, b)
    if ~isa(gc, 'gridcurves')
      gc = plus(-b, gc);
      return
    end
    gc.scalaronly(b)
    for k = 1:numel(gc.curves_)
      gc.curves_{k} = gc.curves_{k} - b;
    end
  end
  
  function gc = mtimes(gc, b)
    if ~isa(gc, 'gridcurves')
      [gc, b] = deal(b, gc);
    end
    gc.scalaronly(b)
    for k = 1:numel(gc.curves_)
      gc.curves_{k} = b*gc.curves_{k};
    end
  end
  
  function gc = mrdivide(gc, b)
    gc = rdivide(gc, b);
  end
  
  function n = numel(gc, varargin)
    % Overloaded for subsref.
    
    n = numel(gc.curves_, varargin{:});
  end
    
  function gc = rdivide(gc, b)
    if isa(gc, 'gridcurves')
      gc.scalaronly(b)
      gc = mtimes(gc, 1/b);
    else
      [gc, b] = deal(b, gc);
      gc.scalaronly(b)
      for k = 1:numel(gc.curves_)
        gc.curves_{k} = b./gc.curves_{k};
      end
    end
  end

  function out = plot(gc, varargin)
    washold = ishold;
    ah = newplot;
    
    gctag = sprintf('gridcurve_%s', num2hex(rand));
    hold on
    
    [gargs, pargs] = cmtplot.gridargs(varargin{:});
    for k = 1:numel(gc.curves_)
      zg = gc.curves_{k};
      line(real(zg), imag(zg), pargs{:}, gargs{:}, 'tag', gctag)
    end
    
    if ~washold
      hold off
    end
    
    if nargout
      out = findobj(ah, 'tag', gctag);
    end
  end
  
  function gc = plus(gc, b)
    if ~isa(gc, 'gridcurves')
      [gc, b] = deal(b, gc);
    end
    gc.scalaronly(b)
    for k = 1:numel(gc.curves_)
      gc.curves_{k} = gc.curves_{k} + b;
    end
  end
    
  function varargout = subsref(gc, S)
    % Provide C(j) or C(j,k) access to curve cell array.
    % Why? See gridcurves help.
    
    switch S(1).type
      case {'()', '{}'}
        if S(1).type(1) == '('
          S(1).type = '{}';
        end
        [varargout{1:nargout}] = subsref(gc.curves_, S);
        
      otherwise
        [varargout{1:nargout}] = builtin('subsref', gc, S);
    end
  end
  
  function gc = uminus(gc)
    for k = 1:numel(gc.curves_)
      gc.curves_{k} = -gc.curves_{k};
    end
  end
end

methods(Access=private)
  function scalaronly(~, b)
    if ~isa(b, 'double') || numel(b) ~= 1
      error('CMT:NotDefined', 'Operation only defined for scalar values.')
    end    
  end
end

end
