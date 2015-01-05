classdef (Abstract) zgrid
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
    region
end

properties (SetAccess = protected)
    dataSource = {}
    dataImage = {}
end

methods
  function g = zgrid(r)
    if ~nargin
      return
    end
    
    g.region = r;
  end
  
  function disp(g)
      fprintf('grid object (%s)\n',class(g))
  end
  
    
end

methods (Abstract)
    w = apply(g,f)    % apply function f to grid g    
    out = plot(g,varargin)   % plot the grid
end
  
end


