classdef cmtplot < cmtobject
% CMTPLOT collects common CMT plot tasks.
% 
% This is probably more convoluted than it needs to be.

% This file is a part of the CMToolbox.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.
% Written by Everett Kropf, 2014.

properties(Constant)
  % colors
  black = 'k'
  grey = 0.75*[1, 1, 1]
  none = 'none'
  
  % line things
  defaultlinewidth = 0.5;
  smoothingon = {'linesmoothing', 'on'}
end

methods
    function p = cmtplot
        % This is here to allow get/set functionality.
        
        get(p, plotset);
    end
end

methods(Static)
  function args = closedcurveargs()
    args = [{'color', cmtplot.cccolor, ...
             'linewidth', cmtplot.cclinewidth}, ...
             cmtplot.smoothingon];
  end
  
  function c = cccolor()
    c = cmtplot.black;
  end
  
  function c = cclinewidth()
    c = cmtplot.defaultlinewidth;
  end
  
  function args = fillargs()
    args = {cmtplot.fillcolor, 'edgecolor', cmtplot.filledgecolor};
  end
  
  function c = fillcolor()
    c = cmtplot.grey;
  end
  
  function c = filledgecolor()
    c = cmtplot.none;
  end
  
  function args = gridargs()
    args = [{'color', cmtplot.gridcolor, ...
             'linewidth', cmtplot.gridlinewidth}, ...
             cmtplot.smoothingon];
  end
  
  function c = gridcolor()
    c = cmtplot.grey;
  end
  
  function c = gridlinewidth()
    c = cmtplot.defaultlinewidth;
  end
  
  function [args, gargs] = pullgridargs(arglist)
    % Separate grid args from cell array arglist. Must have an even number of
    % entries of the {'name', value} pair form.
    
    n = numel(arglist);
    if mod(n, 2)
      error('CMT:InvalidArgument', 'Expected name/value pairs.')
    end
    
    recognized = {'nrad', 'ncirc'};
    gargs = {[], []};
    idx = 1:n;
    for k = 1:2:n-1
      try
        match = validatestring(arglist{k}, recognized);
      catch err
        if strcmp(err.identifier, ...
                  'MATLAB:unrecognizedStringChoice')
          continue
        end
        rethrow(err)
      end
      
      switch match
        case 'nrad'
          gargs{1} = arglist{k+1};
        case 'ncirc'
          gargs{2} = arglist{k+1};
        otherwise
          error('CMT:BadThings', 'This should never happen.')
      end
      
      idx = idx(idx ~= k & idx ~= k+1);
    end
    if ~isempty(idx)
      args = arglist{idx};
    else
      args = {};
    end
  end
  
  function whitefigure(fig)
    if ~nargin
      fig = gcf;
    end
    set(fig, 'color', 'white');
  end
end

end
