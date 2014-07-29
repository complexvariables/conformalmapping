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
    numCollPts              % Number of collcation points.
    kernSolMethod           % Solver method.
    newtonTol               % Newton iteration tolerance.
    trace                   % Print out solution trace information.
    numFourierPts           % Default size of FFT to employ.
end

properties(Access=protected)
    proplist = ...
        {
        'numCollPts', 512, [], '[ integer {512} ]'
        'kernSolMethod', 'auto', [], '[ backslash | orth_resid | {auto} ]'
        'newtonTol', 10*eps(2*pi), [], '[ scalar double {10*eps(2*pi)} ]'
        'trace', false, [], '[ true | {false} ]'
        'numFourierPts', 256, [], '[ integer {256} ]'
        }
end

methods
    function opt = szset(varargin)
        opt = opt@optset(varargin{:});
    end
end
    
end
