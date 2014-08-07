function R = interior(C)
% INTERIOR creates a bounded region with boundary C.
%
% R = exterior(C)
% Creates interior region R bounded by closed curve C.
%
% See also closedcurve, region.

% This file is a part of the CMToolbox.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.

if ~isa(C, 'closedcurve')
    error('CMT:InvalidArgument', ...
        'Function argument must be a closed curve.')
end

R = region(C);
