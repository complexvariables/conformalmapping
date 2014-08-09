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
    lineSmoothing
end

properties(Access=protected)
    proplist = { ...
        'lineWidth', 0.5, @isnumeric, '[ double {0.5} ]'
        'lineSmoothing', 'on', ...
            @plotset.isOnOff, '[ {on} | off ]'
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
