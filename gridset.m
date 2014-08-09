classdef gridset < optset
%GRIDSET holds grid preferences.

% This file is a part of the CMToolkit.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.

properties
    gridType
    numRadialLines
    numCircularLines
    numLevels
end

properties(Access=protected)
    proplist = { ...
        'gridType', 'polar', [], '[ string {polar} | carleson ]'
        'numRadialLines', 20, [], '[ positive integer {20} ]'
        'numCircularLines', 5, [], '[ positive integer {5} ]'
        'numLevels', 5, [], '[ positive integer {5} ]'
    };
end

methods
    function opts = gridset(varargin)
        opts = opts@optset(varargin{:});
    end
    
    function opts = set.gridType(opts, value)
        opts.gridType = lower(value);
    end
end

end
