classdef plotset < optset
% PLOTSET is plot settings class.
%
% See also optset.

% This file is a part of the CMToolkit.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.
% Written by Everett Kropf, 2014.

properties
    lineWidth
    lineColor
    lineSmoothing
    gridColor
    regionColor
end

properties(Access=protected)
    proplist = { ...
        'lineWidth', 2, @isnumeric, '[ double {0.5} ]'
        'lineColor', cmtplot.ver84(1,:), [], '[ valid colorspec ]'
        'lineSmoothing', 'on', ...
            @plotset.isOnOff, '[ {on} | off ]'
        'regionColor',cmtplot.ver84(3,:), [], '[ valid colorspec ]'
        'gridColor', cmtplot.ver84([2 4],:), [], '[ 2 x valid colorspec ]'
    }
end

methods
    function opt = plotset(varargin)
        opt = opt@optset(varargin{:});
    end
    
    function opt = set.lineSmoothing(opt, value)
        opt.lineSmoothing = lower(value);
    end
end

methods(Static, Hidden)
    function tf = isOnOff(s)
        tf = any(strcmpi(s, {'on', 'off'}));
    end
end

end
