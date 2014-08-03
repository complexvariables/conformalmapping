classdef szset < optset
% SZSET is the options structure for the szego kernel.
%
% opts = szset('name', value, ...)
%   Creates option structure for szego via name/value pairs.
%
% Properties:
%   numCollPts        -- Number of collcation points.
%   kernSolMethod     -- Solver method.
%   newtonTol         -- Newton iteration tolerance.
%   trace             -- Print out solution trace information.
%   numFourierPts     -- Default size of FFT to employ.
%
% Methods:
% defaults(szset)
%   Shows properties which may be set along with defaults.
%
% See also szego.

% This file is a part of the CMToolkit.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.
% Written by Everett Kropf, 2014.

properties
    confCenter              % Conformal center.
    numCollPts              % Number of collcation points.
    kernSolMethod           % Solver method.
    newtonTol               % Newton iteration tolerance.
    trace                   % Print out solution trace information.
    numFourierPts           % Default size of FFT to employ.
end

properties(Access=protected)
    proplist = ...
        {
        'confCenter', 0, ...
            @szset.isScalarDouble, '[ scalar double {0} ]'
        'numCollPts', 512, ...
            @szset.isPosInteger, '[ integer {512} ]'
        'kernSolMethod', 'auto', ...
            @szset.isOkMethod, '[ backslash | orth_resid | {auto} ]'
        'newtonTol', 10*eps(2*pi), ...
            @szset.isScalarDouble, '[ scalar double {10*eps(2*pi)} ]'
        'trace', false, [], '[ true | {false} ]'
        'numFourierPts', 256, ...
            @szset.isPosInteger, '[ integer {256} ]'
        }
end

methods
    function opt = szset(varargin)
        opt = opt@optset(varargin{:});
    end
end

methods(Static, Hidden)
    function tf = isPosInteger(x)
        tf = (x - fix(x)) == 0 & x > 0;
    end
    
    function tf = isOkMethod(s)
        tf = any(strcmp(s, {'backslash', 'orth_resid', 'auto'}));
    end
    
    function tf = isScalarDouble(x)
        tf = numel(x) == 1 & isa(x, 'double');
    end
end
    
end
