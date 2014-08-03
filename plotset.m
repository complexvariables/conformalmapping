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
    linewidth
    linesmoothing
end

properties(Access=protected)
    proplist = { ...
        'linewidth', 0.5, @isnumeric, '[ double {0.5} ]'
        'linesmoothing', 'on', ...
            @plotset.isOnOff, '[ {on} | off ]'
    }
end

methods
    function opt = plotset(varargin)
        opt = opt@optset(varargin{:});
    end
end

methods(Static, Hidden)
    function tf = isOnOff(s)
        tf = any(strcmpi(s, {'on', 'off'}));
    end
end

end
