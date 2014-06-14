classdef szset < optset
% SZSET is the options structure for the szego kernel.
%
% opts = szset('name', value, ...)
%   Creates option structure for szego via name/value pairs.
%
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
  nS                              % Number of collcation points.
  method                          % Solver method.
  trace                           % Print out solution trace information.
  nF                              % Default size of FFT to employ.
end

properties(Access=protected)
  proplist = ...
    {
      'nS', 128, [], '[ integer {128} ]'
      'method', 'auto', [], '[ backslash | orth_resid | {auto} ]'
      'trace', false, [], '[ true | {false} ]'
      'nF', 512, [], '[ integer {512} ]'
    }
end

methods
  function opt = szset(varargin)
    opt = opt@optset(varargin{:});
  end
end

end
