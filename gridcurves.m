classdef gridcurves
% GRIDCURVES class holds grid curves for regions.
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
  
  function disp(gd)
    fprintf('gridcurves object:\n\n')
    fprintf('  with %d gridlines.\n\n', numel(gd))
  end
  
  function n = numel(gc, varargin)
    n = numel(gc.curves_, varargin{:});
  end
  
  function out = plot(gc, varargin)
    washold = ishold;
    
    gctag = sprintf('gridcurve_%s', num2hex(rand));
    hold on
    for k = 1:numel(gc.curves_)
      zg = gc.curves_{k};
      args = plotdef.gridargs;
      plot(real(zg), imag(zg), args{:}, varargin{:}, 'tag', gctag)
    end
    
    if ~washold
      hold off
    end
    
    if nargout
      out = findobj(gca, 'tag', gctag);
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
end

end
